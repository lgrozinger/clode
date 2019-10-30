(defsystem "clode"
  :version "0.1.0"
  :author "Lewis Grozinger"
  :license ""
  :depends-on ("gsll")
  :components ((:module "src"
                :components
                ((:file "main"))))
  :description "Interface to the GSL ODE solver."
  :in-order-to ((test-op (test-op "clode/tests"))))

(defsystem "clode/tests"
  :author "Lewis Grozinger"
  :license ""
  :depends-on ("clode"
               "rove")
  :components ((:module "tests"
                :components
                ((:file "main"))))
  :description "Test system for clode"
  :perform (test-op (op c) (symbol-call :rove :run c)))
