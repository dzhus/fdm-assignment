#lang scheme

(require srfi/1
         srfi/43
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
;; initial is a vector of initial values in nodes
(define (solve-heat conductivity
                    initial
                    left-value right-value
                    time-step space-step
                    [layer-length 1])
  ;; We need equations for k from 2 to K-2
  (let ((space-intervals (inexact->exact (round (/ layer-length space-step)))))
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
           (v (vector-map (lambda (k x) (+ x (/ (initial (+ k 1)) time-step)))
                          (vector-append (vector (* r left-value))
                                         (make-vector (- eq-count 2) 0)
                                         (vector (* r right-value))))))
      (vector-append (vector left-value)
                     (solve-tridiagonal A v)
                     (vector right-value)))))

(define-struct grid-point (x y bound value))

(define (distance x1 y1 x2 y2)
  (sqrt (+ (sqr (- x1 x2) ) (sqr (- y1 y2)))))

;; Build a matrix of point structures which form a regular grid
;; withing a bounding box with dimensions specified by `box-x` and
;; `box-y`. Points with coordinates which do not satisfy
;; `body-predicate?` are replaced with #f. For the rest points,
;; `bound` and `initial` structure fields are set to the results of
;; `boundary` and `initial` called with point coordinates,
;; respectively.
(define (make-grid box-x box-y hx hy body-predicate? boundary initial)
  (let ((eps (distance hx hy 0 0)))
    (build-matrix (add1 (inexact->exact (floor (/ box-y hy))))
                  (add1 (inexact->exact (floor (/ box-x hx))))
                  (lambda (i j)
                    (let ((x (* j hx))
                          (y (* i hy)))
                      (if (body-predicate? x y eps)
                          (make-grid-point x y
                                           (boundary x y eps)
                                           (initial x y))
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
        ((= y 6) 50)
        (else #f)))

(define (initial x y) 0)

;; Solve 2D heat problem
(define (solve-2d-heat grid [hx 0.25] [hy 0.25] [dt 25] [conductivity 1])
  (dump-grid-to-file grid "out-initial.txt")
  ;; X
  (for-each (lambda (i)
              (let* ((row (build-vector (matrix-cols grid)
                                        (lambda (j) (matrix-ref grid i j))))
                     (left-i (vector-skip not row))
                     (right-i (vector-skip-right not row))
                     (left-point (vector-ref row left-i))
                     (right-point (vector-ref row right-i)))
                (let ((solution
                       (solve-heat (* 2 conductivity)
                                   (lambda (k) (grid-point-value
                                           (matrix-ref grid i (+ k left-i))))
                                   (grid-point-bound left-point)
                                   (grid-point-bound right-point)
                                   dt hx
                                   (- (grid-point-x right-point)
                                      (grid-point-x left-point)))))
                  (vector-for-each
                   (lambda (j x) (grid-value-update! grid i (+ j left-i) x))
                   solution))))
            (iota (matrix-rows grid)))
  (dump-grid-to-file grid "out-x.txt")
  ;; Y
  (for-each (lambda (i)
              (let* ((row (build-vector (matrix-rows grid)
                                        (lambda (j) (matrix-ref grid j i))))
                     (left-i (vector-skip not row))
                     (right-i (vector-skip-right not row))
                     (left-point (vector-ref row left-i))
                     (right-point (vector-ref row right-i)))
                (let ((solution
                       (solve-heat (* 2 conductivity)
                                   (lambda (k) (grid-point-value
                                           (matrix-ref grid (+ k left-i) i)))
                                   (grid-point-bound left-point)
                                   (grid-point-bound right-point)
                                   dt hy
                                   (- (grid-point-y right-point)
                                      (grid-point-y left-point)))))
                  (vector-for-each
                   (lambda (j x) (grid-value-update! grid (+ j left-i) i x))
                   solution))))
            (iota (matrix-cols grid)))
  (dump-grid-to-file grid "out.txt")
  grid)

(define (run grid hx hy time-step until-time)
  (define (iteration initial-grid at-time)
    (if (< at-time until-time)
        (iteration (solve-2d-heat initial-grid hx hy time-step) (+ at-time time-step))
        at-time))
  (iteration grid 0))

(let ([hx 0.05]
      [hy 1]
      [dt 0.5]
      [until 5])
  (run (make-grid 8 6 hx hy in-body? boundary initial)
       hx hy dt until))
