(defpackage clode/tests/main
  (:use :cl
        :clode
        :rove))
(in-package :clode/tests/main)

;; NOTE: To run this test file, execute `(asdf:test-system :clode)' in your Lisp.

(deftest test-target-1
  (testing "should (= 1 1) to be true"
    (ok (= 1 1))))
