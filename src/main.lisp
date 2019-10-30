(defpackage clode
  (:use :cl :cffi)
  (:import-from :alexandria :flatten :with-gensyms)
  (:export :*lotka-volterra*))
(in-package :clode)

(define-foreign-library libgsl
  (t (:default "libgslcblas")))

(use-foreign-library libgsl)

;; System definition as follows:
;; ((variable . update-form) | (parameter . value))

;; E.g. Lotka-Volterra
(defvar *lotka-volterra* '(x (- (* ɑ x) (* β x y))
                           y (- (* β x y) (* ɣ y))
                           ɑ 10.0
                           β 0.1
                           ɣ 10.0)
  "An example ODE system for Lotka-Volterra dynamics.")

(defun variables (in-system)
  "Collect the variables (including parameters) for the system."
  (loop for (var form) on in-system by #'cddr collect var))

(defun dependent-variables (in-system)
  "Collect the dependent variables for the system."
  (loop for variable in (variables in-system)
     when (not (atom (equation variable in-system)))
     collect variable))

(defun parameters (in-system)
  "Collect the parameters for the system."
  (loop for variable in (variables in-system)
     when (atom (equation variable in-system))
     collect variable))

(defun equation (for-var in-system)
  "Get the equation for the delta of a variable."
  (getf in-system for-var))

(defun equations (in-system)
  "Get all the equations for the system."
  (loop for variable in (variables in-system)
       collect (equation variable in-system)))


;; low level CFFI for libgsl ordinary differential equations
(defcstruct gsl-odeiv2-system
  (function :pointer)
  (jacobian :pointer)
  (dimension :uint)
  (params :pointer))

;; GSL driver generator function
(defcfun "gsl_odeiv2_driver_alloc_y_new" :pointer
  (sys (:pointer (:struct gsl-odeiv2-system)))
  (step-type :pointer)
  (hstart :double)
  (epsabs :double)
  (epsrel :double))

(defun symbol-to-mem-aref (form symbol pointer index)
  "Replace all occurrences of symbol in form with (mem-aref pointer :double index)."
  (cond ((and (symbolp form) (eql form symbol))
         (list 'mem-aref pointer :double index))
        ((listp form)
         (loop for sub-form in form
            collect (symbol-to-mem-aref sub-form symbol pointer index)))
        (t form)))

(defun system-to-mem-aref (system y-pointer p-pointer)
  "Replace occurrences of symbols with appropriate mem-aref forms."
  (let ((result system)
        (ys (dependent-variables system))
        (ps (parameters system)))
    (loop for y in ys for i from 0 below (length ys)
       do (setf result (loop for form in result
                          collect (symbol-to-mem-aref form y y-pointer i))))
    (loop for p in ps for i from 0 below (length ps)
       do (setf result (loop for form in result
                       collect (symbol-to-mem-aref form p p-pointer i))))

    result))

(defun nested-replace (new old sequence)
  (cond ((eql sequence old) new)
        ((listp sequence) (loop for sub in sequence collect (nested-replace new old sub)))
        (t sequence)))

(defmacro defode ((name) &body system-definition)
  "Construct a C CFFI callback to be used for GSL ODE methods that implements the system definition."
  (with-gensyms (time y dydt params)
    (let ((transformed-system (system-to-mem-aref system-definition y params)))
      `(defcallback ,name :int ((,time :double)
                                (,y (:pointer :double))
                                (,dydt (:pointer :double))
                                (,params :pointer))
         (declare (ignorable ,time ,params))
         ,@(loop for (variable update) on transformed-system by #'cddr
              when (member y variable)
              collect (list 'setf (nested-replace dydt y variable) update))

         0))))

(defmacro with-ode-system ((name &rest system) &body body)
  "Construct a gsl_odeiv2_system for system, bind to name and execute body."
  (with-gensyms (foreign-params odefun)
    `(progn
       (defode (,odefun) ,@system)
       (with-foreign-objects ((,name '(:struct gsl-odeiv2-system)) (,foreign-params :double ,(length (parameters system))))
         ;; populate the foreign allocated  parameter array
         ,@(loop for i from 0 below (length (parameters system)) for p in (parameters system)
              collect (list 'setf (list 'mem-aref foreign-params :double i) (coerce (getf system p) 'double-float)))
         
         ;; set the ODE callback and Jacobian for GSL ODE system foreign struct              
         (setf (foreign-slot-value ,name '(:struct gsl-odeiv2-system) 'function) (callback ,odefun)
               (foreign-slot-value ,name '(:struct gsl-odeiv2-system) 'jacobian) (null-pointer)
               (foreign-slot-value ,name '(:struct gsl-odeiv2-system) 'dimension) ,(length (dependent-variables system))
               (foreign-slot-value ,name '(:struct gsl-odeiv2-system) 'params) ,foreign-params)

         ,@body))))       


         

        
  

