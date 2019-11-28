(defpackage clode
  (:use :cl :cffi)
  (:import-from :alexandria :flatten :with-gensyms)
  (:export :*lotka-volterra*
           :*rk23* :*rk4* :*rk45* :*rkck* :*rk8pd* :*msadams*
           :with-ode-system
           :integrate
           :with))
(in-package :clode)

;; put your libgsl shared library here
(define-foreign-library libgsl
  (t (:default "libgsl")))

(use-foreign-library libgsl)

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

;; GSL driver set minimum step
(defcfun "gsl_odeiv2_driver_set_hmin" :int
  (d :pointer)
  (hmin :double))

;; GSL driver free function
(defcfun "gsl_odeiv2_driver_free" :int
  (d :pointer))

(defun withp (token)
  "Is the token 'with'"
  (and
   (symbolp token)
   (string= (symbol-name token) (symbol-name 'with))))

(defun fromp (token)
  "Is the token 'from'"
  (and
   (symbolp token)
   (string= (symbol-name token) (symbol-name 'from))))

(defun byp (token)
  "Is the token 'by'"
  (and
   (symbolp token)
   (string= (symbol-name token) (symbol-name 'by))))

       
(defun equations-clause (definition)
  (loop
     for token in definition
     until (or (withp token) (fromp token))
     collect token))

(defun variables-clause (definition)
  (loop
     with flag = NIL
     for token in definition
     until (fromp token)
     when flag
     collect token
     do (setf flag (or flag (withp token)))))

(defun from-to-clause (definition)
  (loop
     with flag = NIL
     for token in definition
     until (byp token)
     when flag
     collect token
     do (setf flag (or flag (fromp token)))))

(defun by-clause (definition)
  (loop
     with flag = NIL
     for token in definition
     when flag
     collect token
     do (setf flag (or flag (byp token)))))

(defun recursive-replace (new old form)
  "Replace (recursively) all instances of old with new in form."
  (flet ((replacer (x)
           (cond ((eql x old) new)
                 ((listp x) (recursive-replace new old x))
                 (t x))))
    (if (listp form)
        (mapcar #'replacer form)
        (replacer form))))

(defun state-variable-index (definition)
  (loop
     for (variable equal-sign form) on (equations-clause definition) by #'cdddr
     for i = 0 then (+ 1 i)
     collect variable
     collect i))

(defun parameter-index (definition)
  (loop
     for (parameter equal-sign form) on (variables-clause definition) by #'cdddr
     with i = 0
     with index = '()
     do (unless (member parameter (state-variable-index definition))
          (setf index (append index (list parameter i)))
          (setf i (+ 1 i)))
     finally (return index)))

(defun value-of-variable (variable definition)
  (let ((clause (variables-clause definition)))
    (nth (+ 2 (position variable clause)) clause)))

(defun substitute-symbol-for-foreign (variables foreign-pointer equations)
  (loop
     for (variable index) on variables by #'cddr
     do (setf equations
              (loop
                 for (v _ form) on equations by #'cdddr
                 collect v
                 collect '=
                 collect (recursive-replace
                          `(mem-aref ,foreign-pointer :double ,index)
                          variable
                          form)))
     finally (return equations)))

(defmacro integrate (&rest definition)    
  (with-gensyms (callback
                 system-name
                 f-time
                 f-y
                 f-dydt
                 f-params)
    (let* ((ps (parameter-index definition))
           (ys (state-variable-index definition))
           (equations
            (substitute-symbol-for-foreign
              ps
              f-params
              (substitute-symbol-for-foreign
               ys
               f-y
               (equations-clause definition)))))

      `(progn
         (defcallback ,callback :int ((,f-time :double)
                                      (,f-y (:pointer :double))
                                      (,f-dydt (:pointer :double))
                                      (,f-params (:pointer :double)))
           (declare (ignorable ,f-time ,f-params))
           ,@(loop
                for (lhs _ rhs) on equations by #'cdddr
                collect `(setf (mem-aref ,f-dydt :double
                                         ,(getf ys lhs))
                               ,rhs))
           0)

         (with-foreign-objects ((,system-name '(:struct gsl-odeiv2-system))
                                 (,f-params :double ,(/ (length ps) 2)))
           ;; populate the foreign allocated parameter array
           ,@(loop
                for (parameter i) on ps by #'cddr
                collect `(setf (mem-aref ,f-params :double ,i)
                               ,(value-of-variable parameter definition)))

           ;; set the ODE callback and Jacobian for GSL ODE system foreign struct
           (with-foreign-slots ((function jacobian dimension params)
                                ,system-name (:struct gsl-odeiv2-system))
           
             (setf function (callback ,callback)
                   jacobian (null-pointer)
                   dimension ,(/ (length ys) 2)
                   params ,f-params)
            
             (integrate-1 ,system-name
                          '(,@(loop
                                 for (variable _) on ys by #'cddr
                                 collect (value-of-variable variable definition)))
                          :start ,(first (from-to-clause definition))
                          :end ,(third (from-to-clause definition))
                          :steps ,(first (by-clause definition)))))))))

(defun integrate-1 (gsl-system initial-state
                  &key
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
               (if (not (eq status 0))
                   (format t "~&The integrator failed with status ~A~%" status))
             collect
               (cons (mem-ref time :double)
                     (loop for i below d collect (mem-aref y :double i)))
             finally
               (gsl-odeiv2-driver-free driver)))))))
             
;; E.g. Lotka-Volterra
(defun example ()
  "An example ODE system for Lotka-Volterra dynamics."
  (integrate x = (- (* ɑ x) (* β x y))
             y = (- (* β x y) (* ɣ y))
             with x = 1.d2 y = 1.d1 ɑ = 1.d1 β = 1.d-1 ɣ = 1.d1
             from 0.d0 to 1.d2 by 10000))
      
      
  
  
                                        


