(ns kyleerhabor.elliptic-curve-discrete-logarithm-problem-collision-algorithm.core
  (:require
   [clojure.math :refer [ceil sqrt]]
   [clojure.tools.logging :as log]))

(defn point-name [{:keys [x y]
                   :as p}]
  (if (nil? p)
    "O"
    (str "(" x ", " y ")")))

(defn theorem-name [{:keys [a b c p pp pq]}]
  (str "C: y^2 = x^3 + " a "x^2 + " b "x + " c ", p = " p ", P = " (point-name pp) ", Q = " (point-name pq)))

(defn square [x]
  (* x x))

(defn mod-inverse [a p]
  (long (.modInverse (biginteger a) (biginteger p))))

(defn add-point [{:keys [a p]}
                 lambda
                 {px :x
                  py :y}
                 {qx :x}]
  (let [;; ν = y_1 − λx_1 = y_2 − λx_2
        nu (mod (- py (* lambda px)) p)
        ;; x_3 = λ^2 - a - x_1 - x_2
        x (mod (- (square lambda) a px qx) p)
        ;; y_3 = -λx_3 - ν
        y (mod (- (- (* lambda x)) nu) p)]
    {:x x
     :y y}))

(defn add-duplicate [{:keys [a b p]
                      :as ecdlp}
                     {px :x
                      py :y
                      :as pp}
                     pq]
  (let [;; y^2 = f(x) = x^3 + ax^2 + bx + c
        ;; d/dx f(x) = 3x^2 + 2ax + b
        ;; λ = d/dx y | P_0
        ;;   = f'(x_1) / 2y_1
        ;;   = d/dx f(x) / 2y_1
        ;;   = (3x_1^2 + 2ax_1 + b) / 2y_1
        lambda (mod (*
                      (+ (* 3 (square px)) (* 2 a px) b)
                      (mod-inverse (* 2 py) p)) p)]
    (add-point ecdlp lambda pp pq)))

(defn- add-general [{:keys [p]
                     :as ecdlp}
                    {px :x
                     py :y
                     :as pp}
                    {qx :x
                     qy :y
                     :as pq}]
  (let [;; λ = (y_2 - y_1) / (x_2 - x_1)
        lambda (mod (*
                      (- qy py)
                      (mod-inverse (- qx px) p)) p)]
    (add-point ecdlp lambda pp pq)))

(defn add [ecdlp k pp pq]
  (let [point (cond
                (nil? pp) pq
                (nil? pq) pp ;; I think this is dead code for our problem, but I'll leave it in for completeness.
                (= pp pq) (add-duplicate ecdlp pp pq)
                :else (add-general ecdlp pp pq))]
    (log/info (str (theorem-name ecdlp) ": k = " k ", " (point-name pp) " + " (point-name pq) " = " (point-name point)))
    point))

(defn negate [{:keys [x y]}]
  {:x x
   :y (- y)})

(defn theorem [{:keys [p pp pq]
                :as ecdlp}]
  (let [;; 4.28 (Shank's Babystep-Giantstep Algorithm).
        ;;
        ;; Step 1: Let n = ⌈√N⌉ be the smallest integer that is greater than √N.
        ;;
        ;; n doesn't have to equal ⌈√N⌉ for the algorithm to work, but it helps in making it efficient. We can use √p as
        ;; an approximate of √N.
        n (ceil (sqrt p))
        _ (log/info (str (theorem-name ecdlp) ": n = " n))
        baby-step-base (add ecdlp 0 pp nil)
        {:keys [steps latest]} (loop [ks (range 2 (inc n))
                                      r baby-step-base
                                      {:keys [steps]
                                       :as data} {:steps {}
                                                  :latest nil}]
                                 (if-let [k (first ks)]
                                   (let [r (add ecdlp k baby-step-base r)]
                                     (recur (rest ks) r {:steps (assoc steps r k)
                                                         :latest r}))
                                   data))
        giant-step-base (negate latest)
        m (loop [ks (range n)
                 r pq]
            (let [k (first ks)]
              (if-let [a (get steps r)]
                (do
                  (log/info (str (theorem-name ecdlp) ": " a "P = Q - " n " * " k "P"))
                  (long (+ a (* n k))))
                (recur (rest ks) (add ecdlp k r giant-step-base)))))]
    (log/info (str (theorem-name ecdlp) ": " m "P = Q"))
    m))

(defn -main []
  ;; C: y^2 = x^3 + x^2 + x + 1, p = 97, P = (7, 20), Q = (17, 46).
  (assert (= 47 (theorem {:a 1 :b 1 :c 1
                          :p 97
                          :pp {:x 7 :y 20}
                          :pq {:x 17 :y 46}})))

  ;; (a) C: y^2 = x^3 + x^2 + x + 3, p = 103, P = (7, 14), Q = (8, 22).
  (assert (= 76 (theorem {:a 1 :b 1 :c 3
                          :p 103
                          :pp {:x 7 :y 14}
                          :pq {:x 8 :y 22}})))

  ;; (b) C: y^2 = x^3 - 2x^2 + 5x + 6, p = 149, P = (11, 16), Q = (110, 46).
  (assert (= 87 (theorem {:a -2 :b 5 :c 6
                          :p 149
                          :pp {:x 11 :y 16}
                          :pq {:x 110 :y 46}})))
  ;; (c) C: y^2 = x^3 + x^2 + x + 2, p = 10037, P = (8, 7358), Q = (2057, 5437).
  (assert (= 1277 (theorem {:a 1 :b 1 :c 2
                            :p 10037
                            :pp {:x 8 :y 7358}
                            :pq {:x 2057 :y 5437}}))))
