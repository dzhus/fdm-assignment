#lang scheme

(require srfi/43
         srfi/48
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
           (v (vector-map (lambda (k x) (+ x (/ (initial (+ k 2)) time-step)))
                          (vector-append (vector left-value)
                                         (make-vector (- eq-count 2) 0)
                                         (vector right-value)))))
      (vector-append (vector left-value)
                     (solve-tridiagonal A v)
                     (vector right-value)))))

;; Generate a list of solutions for heat equation on every time layer.
;; Return a list of pairs where car is time and cdr is a vector of
;; temperature values along the interval.
(define (solve-heat-problem conductivity
                            initial
                            left-bound right-bound
                            time-step space-step
                            until-time)
  (define (function-to-initial f)
    (lambda (k) (f (* (sub1 k) space-step))))
  (define (iteration initial at-time [acc '#()])
    (if (< at-time until-time)
        (let* ((solution (solve-heat conductivity initial
                                     (left-bound at-time)
                                     (right-bound at-time)
                                     time-step space-step))
               (next-initial
                (lambda (k) (vector-ref solution (sub1 k)))))
          (iteration next-initial (+ at-time time-step)
                     (vector-append acc (vector (cons at-time solution)))))
        acc))
  (iteration (function-to-initial initial) time-step))

(define (simple-print solutions space-step)
  (vector-for-each
   (lambda (k layer)
     (let ((at-time (car layer))
           (solution (cdr layer)))
       (vector-for-each
        (lambda (k value)
          (display (format "~0,2F ~0,2F ~a~%" at-time (* space-step k) value)))
        solution))
     (newline))
   solutions))

(define (meshviewer-print solutions space-step)
  (let ((time-layers (vector-length solutions))
        (space-layers (add1 (inexact->exact (/ space-step)))))
    (define (node-index space-index time-index)
      (+ (add1 space-index) (* time-index space-layers)))
    ;; Return index for element with (node-index s-i t-i) in upper left corner
    (define (element-index space-index time-index)
      (+ (add1 space-index) (* time-index (sub1 space-layers))))

    ;; Print solutions for every time layer, sequentially enumerating
    ;; all nodes
    (define (print-nodes solutions)
      (vector-for-each
       (lambda (time-index layer)
         (let ((at-time (car layer))
               (values (cdr layer)))
           (vector-for-each
            (lambda (space-index value)
              (display (format "~d ~0,2F ~0,2F ~a~%"
                               (node-index space-index time-index)
                               at-time
                               (* space-step space-index)
                               value)))
            values)))
       solutions))
    
    (define (print-elements solutions)
      (define (last-space-index? k)
        (= k (sub1 space-layers)))
      (define (last-time-index? k)
        (= k (sub1 time-layers)))
      (vector-for-each
       (lambda (time-index layer)
         (when (not (last-time-index? time-index))
           (let ((at-time (car layer))
                 (values (cdr layer)))
             (vector-for-each
              (lambda (space-index value)
                (when (not (last-space-index? space-index))
                  (let* ((n1 (node-index space-index time-index))
                         (n2 (add1 n1))
                         (n3 (node-index space-index (add1 time-index)))
                         (n4 (add1 n3)))
                    (display (format "~d ~d ~d ~d ~d~%"
                                     (element-index space-index time-index)
                                     n1 n2 n3 n4)))))
                values))))
       solutions))
    
    (display "# Nodes\n")
    (let ((node-count (* time-layers space-layers)))
      (display (format "~a ~a ~a ~a~%" node-count 2 1 "u")))
    (print-nodes solutions)
    (display "# Elements\n")
    (let ((element-count (* (sub1 time-layers) (sub1 space-layers))))
      (display (format "~a ~a ~a~%" element-count 4 0)))
    (print-elements solutions)))

(meshviewer-print
 (solve-heat-problem
  0.02
  (lambda (x)
    (* 100 (sin (* pi x))))
  (lambda (t) 0)
  (lambda (t) (* 10 (sin t)))
  0.1
  0.1
  (command-line
   #:program "temp"
   #:args (until)
   (string->number until)))
 0.1)
