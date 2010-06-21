#lang scheme

(require srfi/1
         srfi/43
         srfi/48
         ;; Windows:
         ;;simple-matrix/matrix
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
;; initial is a function which return initial value in k-th node,
;; where k is 0-based
(define (solve-heat conductivity
                    initial
                    left-value right-value
                    time-step space-step
                    [layer-length 1]
                    [left-flow #f] [right-flow #f])
  ;; We need equations for k from 2 to K-2
  (let ((space-intervals (inexact->exact (round (/ layer-length space-step)))))
    ;; k is 1-based
    (let* ((eq-count (sub1 space-intervals))
           (r (/ conductivity (sqr space-step)))
           (A (build-matrix eq-count eq-count
                            (lambda (k n)
                              (cond ((= k n)
                                     (+ (/ time-step)
                                        (* (if (or (and left-flow (= k 0))
                                                   (and right-flow (= k (sub1 eq-count))))
                                               1
                                               2) r)))
                                    ((= (abs (- k n)) 1)
                                     (- r))
                                    (else 0)))))
           (v (vector-map (lambda (k x) (+ x (/ (initial (+ k 2)) time-step)))
                          (vector-append (vector (if left-flow
                                                     (* r space-step left-flow)
                                                     (* r left-value)))
                                         (make-vector (- eq-count 2) 0)
                                         (vector (if right-flow
                                                     (* r space-step right-flow)
                                                     (* r right-value)))))))
      (let ((solution (solve-tridiagonal A v)))
      (vector-append (vector (if left-flow
                                 (- (vector-ref solution 0) (* space-step (- left-flow)))
                                 left-value))
                     solution
                     (vector (if right-flow
                                 (- (vector-ref solution (sub1 (vector-length solution))) (* space-step (- right-flow)))
                                 right-value)))))))

(define-struct grid-point (x y bound value flow))

(define (distance x1 y1 x2 y2)
  (sqrt (+ (sqr (- x1 x2) ) (sqr (- y1 y2)))))

;; Build a matrix of point structures which form a regular grid
;; withing a bounding box with dimensions specified by `box-x` and
;; `box-y`. Points with coordinates which do not satisfy
;; `body-predicate?` are replaced with #f. For the rest points,
;; `bound`, `initial` and `flow` structure fields are set to the
;; results of `boundary`, `initial` and `flow` arguments called with
;; point coordinates, respectively. If `bound` result is non-nil, then
;; initial is set to it instead.
(define (make-grid box-x box-y hx hy body-predicate? boundary initial flow)
  (let ((eps (distance hx hy 0 0)))
    (build-matrix (add1 (inexact->exact (floor (/ box-y hy))))
                  (add1 (inexact->exact (floor (/ box-x hx))))
                  (lambda (i j)
                    (let ((x (* j hx))
                          (y (* i hy)))
                      (if (body-predicate? x y eps)
                          (let ((bound (boundary x y eps)))
                            (make-grid-point x y
                                             bound
                                             (if bound bound (initial x y))
                                             (flow x y)))
                            #f))))))

;; Set `value` field of point structure in grid to `new-point-value`
(define (grid-value-update! grid i j new-point-value)
  ;(display (format "~a → ~a" (grid-point-value (matrix-ref grid i j)) new-point-value))
  (matrix-set! grid i j
               (struct-copy grid-point
                            (matrix-ref grid i j)
                            [value new-point-value])))

(define (display-grid grid)
  (for/list ((p (in-matrix grid)))
            (when p
              (display (format "~a ~a ~a\n"
                               (grid-point-x p)
                               (grid-point-y p)
                               (or (grid-point-value p) 0)))))
  (void))

(define (dump-grid-to-file grid filename)
  (with-output-to-file filename
    (lambda ()
      (display-grid grid))
    #:exists 'replace))

(define (in-body? x y [eps 0.1])
  (or (and (<= y 6) (<= x 5))
      (and (<= y 3) (<= x 8))
      (<= (distance x y 5 3) 3)))

(define (boundary x y [eps 0.1])
  (cond ((= x 0) 200)
        ((or (= x 8)
             (and (>= x 5) (>= y 3) (<= (abs (- (distance x y 5 3) 3)) eps))) 100)
        ((= y 0) 50)
        (else #f)))

(define (flow x y)
  (cond ((= y 6) 0)
        ((and (= y 0) (= x 4)) 500)
        (else #f)))

(define (initial x y) 0)

;; Solve 2D heat problem
(define (solve-2d-heat grid hx hy dt [conductivity 1])
    ;; X
  (for-each (lambda (i)
              (let* ((row (build-vector (matrix-cols grid)
                                        (lambda (j) (matrix-ref grid i j))))
                     (left-i (vector-skip not row))
                     (right-i (vector-skip-right not row))
                     (left-point (vector-ref row left-i))
                     (right-point (vector-ref row right-i)))
                (let ((solution
                       (solve-heat conductivity
                                   (lambda (k) (grid-point-value
                                           ;; Initial function is 1-based
                                           (matrix-ref grid i (sub1 (+ k left-i)))))
                                   (grid-point-bound left-point)
                                   (grid-point-bound right-point)
                                   dt hx
                                   (- (grid-point-x right-point)
                                      (grid-point-x left-point))
                                   (grid-point-flow left-point)
                                   (grid-point-flow right-point))))
                  (vector-for-each
                   (lambda (j x) (grid-value-update! grid i (+ j left-i) x))
                   solution))))
            (iota (- (matrix-rows grid) 2) 1))
  ;; Y
  (for-each (lambda (i)
              (let* ((row (build-vector (matrix-rows grid)
                                        (lambda (j) (matrix-ref grid j i))))
                     (left-i (vector-skip not row))
                     (right-i (vector-skip-right not row))
                     (left-point (vector-ref row left-i))
                     (right-point (vector-ref row right-i)))
                (let ((solution
                       (solve-heat conductivity
                                   (lambda (k) (grid-point-value
                                           (matrix-ref grid (sub1 (+ k left-i)) i)))
                                   (grid-point-bound left-point)
                                   (grid-point-bound right-point)
                                   dt hy
                                   (- (grid-point-y right-point)
                                      (grid-point-y left-point))
                                   (grid-point-flow left-point)
                                   (grid-point-flow right-point))))
                  (vector-for-each
                   (lambda (j x) (grid-value-update! grid (+ j left-i) i x))
                   solution))))
            (iota (- (matrix-cols grid) 2) 1))
  grid)

(define (run grid hx hy time-step until-time dump-file-name)
  (define (solve initial-grid at-time)
    (if (< at-time until-time)
        (solve (solve-2d-heat initial-grid hx hy time-step) (+ at-time time-step))
        initial-grid))
  (dump-grid-to-file (solve grid 0) dump-file-name))

(let ([hx 0.25]
      [hy 0.25]
      [dt 0.05]
      [until 5])
  (run (make-grid 8 6 hx hy in-body? boundary initial flow)
       hx hy dt until
       "out-c.txt"))
