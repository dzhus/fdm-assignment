#lang scheme

(require srfi/43
         (planet wmfarr/simple-matrix:1:0/matrix)
         pyani-lib/matrix
         pyani-lib/linear-eq)

(provide/contract
 [solve-heat
  (number? procedure? number? number? number? number? . -> . vector)])

;; Calculate next layer for solution of heat equation
;;
;; left-value and right-value are boundary values
;;
;; initial must accept integer argument which is 1-based index of node
;; and return the initial value in that node
(define (solve-heat conductivity
                    initial
                    left-value right-value
                    time-step space-step)
  ;; We need equations for k from 2 to K-2
  (let ((space-intervals (inexact->exact (/ space-step))))
    ;; k is 1-based
    (let* ((eq-count (sub1 space-intervals))
           (r (/ conductivity (sqr space-step)))
           (A (build-matrix eq-count eq-count
                            (lambda (k n)
                              (cond ((= k n)
                                     (+ (/ time-step) (* 2 r)))
                                    ((= (abs (- k n)) 1)
                                     (- r))
                                    (else 0)))))
           (v (build-vector eq-count
                            (lambda (k)
                              (let ((k (+ 2 k)))
                                (+ (/ (initial k) time-step)
                                   (cond ((= k 2)
                                          (* r left-value))
                                         ((= k (add1 eq-count))
                                          (* r right-value))
                                         (else 0))))))))
      (vector-append (vector left-value)
                     (solve-tridiagonal A v)
                     (vector right-value)))))

(define (solve-heat-problem conductivity
                            initial
                            left-bound right-bound
                            time-step space-step
                            until-time)
  (define (function-to-initial f)
    (lambda (k) (f (* (sub1 k) space-step))))
  (define (iteration initial at-time acc)
    (if (< at-time until-time)
        (let* ((solution (solve-heat conductivity initial
                                     (left-bound at-time)
                                     (right-bound at-time)
                                     time-step space-step))
               (next-initial
                (lambda (k) (vector-ref solution (sub1 k)))))
          (iteration next-initial (+ at-time time-step)
                     (append acc (list (cons at-time solution)))))
        acc))
  (iteration (function-to-initial initial) time-step '()))

(define (simple-print solutions space-step)
  (for-each
   (lambda (layer)
     (let ((at-time (car layer))
           (solution (cdr layer)))
       (vector-for-each
        (lambda (k value)
          (display (format "~a ~a ~a~n" at-time (* space-step k) value)))
        solution))
     (newline))
   solutions))

(simple-print (solve-heat-problem 0.02
                                  (lambda (x)
                                    (* 100 (sin (* pi x))))
                                  (lambda (t) 0)
                                  (lambda (t) 0)
                                  1
                                  0.1
                                  (command-line
                                   #:program "temp"
                                   #:args (until)
                                   (string->number until)))
              0.1)
