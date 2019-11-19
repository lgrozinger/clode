(defpackage clode
  (:use :cl :cffi)
  (:import-from :alexandria :flatten :with-gensyms)
  (:export :*lotka-volterra*
           :*rk23* :*rk4* :*rk45* :*rkck* :*rk8pd* :*msadams*
           :with-ode-system
           :integrate))
(in-package :clode)

;; put your libgsl shared library here
(define-foreign-library libgsl
  (:unix "/usr/local/lib/libgsl.so.25.0.0")
  (t (:default "libgslcblas")))

(use-foreign-library libgsl)

;; System definition as follows:
;; ((variable . update-form) | (parameter . numberp))

;; E.g. Lotka-Volterra
(defvar *lotka-volterra* '(x (- (* ɑ x) (* β x y))
                           y (- (* β x y) (* ɣ y))
                           ɑ 10.0
                           β 0.1
                           ɣ 10.0)
  "An example ODE system for Lotka-Volterra dynamics.")


;; GSL CFFI DEFINITION

;; GSL stepping methods
(defcvar (*rk23* "gsl_odeiv2_step_rk2") :pointer)
(defcvar (*rk4* "gsl_odeiv2_step_rk4") :pointer)
(defcvar (*rk45* "gsl_odeiv2_step_rkf45") :pointer)
(defcvar (*rkck* "gsl_odeiv2_step_rkck") :pointer)
(defcvar (*rk8pd* "gsl_odeiv2_step_rk8pd") :pointer)
(defcvar (*msadams* "gsl_odeiv2_step_msadams") :pointer)

;; GSL ode system
(defcstruct gsl-odeiv2-system
  (function :pointer)
  (jacobian :pointer)
  (dimension :uint)
  (params :pointer))

(defctype odes (:struct gsl-odeiv2-system))

;; GSL driver generator function
(defcfun "gsl_odeiv2_driver_alloc_y_new" :pointer
  (sys (:pointer odes))
  (step-type :pointer)
  (hstart :double)
  (epsabs :double)
  (epsrel :double))

;; GSL driver application function
(defcfun "gsl_odeiv2_driver_apply" :int
  (d :pointer)
  (time :pointer)
  (next-time :double)
  (y (:pointer :double)))

;; GSL driver free function
(defcfun "gsl_odeiv2_driver_free" :int
  (d :pointer))


(defun variables (in-system)
  "Collect the variables (including parameters) for the system."
  (loop for (var form) on in-system by #'cddr collect var))

(defun dependent-variables (in-system)
  "Collect the dependent variables for the system."
  (loop for variable in (variables in-system)
     when (not (numberp (equation variable in-system)))
     collect variable))

(defun parameters (in-system)
  "Collect the parameters for the system."
  (loop for variable in (variables in-system)
     when (numberp (equation variable in-system))
     collect variable))

(defun equation (for-var in-system)
  "Get the equation for the delta of a variable."
  (getf in-system for-var))

(defun equations (in-system)
  "Get all the equations for the system."
  (loop for variable in (variables in-system)
     collect (equation variable in-system)))

(defun tree-substitute (new old tree)
  "Replace (recursively) all instances of old with new in tree."
  (flet ((replacer (x)
           (cond ((eql x old) new)
                 ((listp x) (tree-substitute new old x))
                 (t x))))
    (mapcar #'replacer tree)))

(defun substitute-variables (variables array-pointer system)
  (let ((result system)
        (dimension (length variables)))
    
    (flet ((array-ref (i) (list 'mem-aref array-pointer :double i)))
      (loop for variable in variables for i below dimension
         do (setf result (tree-substitute (array-ref i) variable result))))

    result))

(defmacro defode ((name) &body system)
  "Define a CFFI callback for the ODE system."
  (with-gensyms (time y dydt params)
    ;; replace dependent-variables with mem-aref forms
    (setf system (substitute-variables (dependent-variables system) y system))

    ;; replace parameters with mem-aref forms
    (setf system (substitute-variables (parameters system) params system))

    ;; replace occurence of the symbol TIME with the gensym'ed symbol
    (setf system (tree-substitute time 'TIME system))

    ;; here is the callback definition
    `(defcallback ,name :int ((,time :double)
                              (,y (:pointer :double))
                              (,dydt (:pointer :double))
                              (,params :pointer))
       (declare (ignorable ,time ,params))
       
       ;; convert the system to a progn of SETF forms
       ,@(loop for (lhs rhs) on system by #'cddr when (member y lhs)
            collect (list 'setf (substitute dydt y lhs) rhs))

       ;; finally return the GSL_SUCCESS code
       0)))

(defmacro with-ode-system ((name &rest system) &body body)
  "Construct a gsl_odeiv2_system for system, bind to name and execute body."
  (with-gensyms (foreign-params odefun)
    (let ((params (parameters system)))
    `(progn
       (defode (,odefun) ,@system)
       (with-foreign-objects ((,name '(:struct gsl-odeiv2-system))
                              (,foreign-params :double ,(length params)))
         ;; populate the foreign allocated parameter array
         ,@(loop for param in params for i below (length params)
              collect
                (list 'setf
                      (list 'mem-aref foreign-params :double i)
                      (coerce (getf system param) 'double-float)))
         
         ;; set the ODE callback and Jacobian for GSL ODE system foreign struct
         (with-foreign-slots ((function jacobian dimension params)
                              ,name (:struct gsl-odeiv2-system))
           
           (setf function (callback ,odefun)
                 jacobian (null-pointer)
                 dimension ,(length (dependent-variables system))
                 params ,foreign-params))
           
         ,@body)))))

(defun integrate (gsl-system initial-state
                  &key
                    (out-file t)
                    (start 0.d0)
                    (end 1.d1)
                    (steps 100)
                    (method *rk8pd*)
                    (abserr 1.d-6)
                    (relerr 1.d-2))

  (let ((d (foreign-slot-value gsl-system 'odes 'dimension)))
    (with-foreign-objects ((y :double d) (time :double))
      (let  ((step-size (/ (- end start) steps))
             (driver (gsl-odeiv2-driver-alloc-y-new
                      gsl-system
                      method
                      1.d-6
                      abserr
                      relerr)))

        ;; set the initial state in the foreign array y
        (loop for i below d
           do (setf (mem-aref y :double i)
                    (coerce (elt initial-state i) 'double-float)))

        ;; set the initial time
        (setf (mem-ref time :double) start)

        ;; termination test
        (flet ((endedp (next-time)
                 (if (> (- end start) 0)
                     (> next-time end)
                     (< next-time end))))
          ;; the integration loop
          (loop
             for next-time = (+ start step-size)
             then (+ (mem-ref time :double) step-size)
             with status = 0
             until (or (endedp next-time) (not (eql status 0)))
             do
               (setf status (gsl-odeiv2-driver-apply driver time next-time y))
               
               (case status
                 (0 (format out-file "~&~,6f ~{~,6f~^ ~}~%"
                            (mem-ref time :double)
                            (loop for i below d collect (mem-aref y :double i))))
                 (otherwise (format t "~&The integrator failed with status ~A~%"
                                    status)
                            (return)))

             finally (gsl-odeiv2-driver-free driver))

          (cons
           (mem-ref time :double)
           (loop for i below d collect (mem-aref y :double i))))))))
             
             
      
      
  
  
                                        


