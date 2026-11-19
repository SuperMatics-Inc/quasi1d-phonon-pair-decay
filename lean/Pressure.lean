import Mathlib.Data.Rat.Defs
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

namespace AdiabaticPressure

theorem covarianceHalf (a da b db : ℚ) (ha : 2*a*da = -1) (hb : 2*b*db = 1) :
    a*da = -1/2 ∧ b*db = 1/2 := by
  constructor <;> linarith

theorem lehmannPositive (wr wl gap : ℚ) (hw : wl ≤ wr) (hg : 0 ≤ gap) :
    0 ≤ (wr-wl)*gap := mul_nonneg (sub_nonneg.mpr hw) hg

theorem lehmannPair (wr wl gap : ℚ) :
    (wr-wl)*gap + (wl-wr)*(-gap) = 2*((wr-wl)*gap) := by ring

theorem nonzeroModeNumerator (u mode omega0 omega1 : ℚ) :
    u*omega1^2*(mode^2+omega0^2) - u*omega0^2*(mode^2+omega1^2)
      = u*mode^2*(omega1^2-omega0^2) := by ring

theorem brownianNormalization (beta u omega : ℚ) :
    beta * (u*omega^2) = u * beta * omega^2 := by ring

theorem factorialStep (j : ℚ) (hj : 1 ≤ j) : 6*j ≤ (2*j)*(2*j+1) := by
  nlinarith [mul_nonneg (show 0 ≤ j by linarith) (show 0 ≤ j-1 by linarith)]

theorem thermalCoefficient (x : ℚ) : (x/2)^2/6 = x^2/24 := by ring

theorem signedFluxError (D c budget : ℚ) (hD : 0 ≤ D) (h1 : -1 ≤ c)
    (h2 : c ≤ 1) (hb : 2*D ≤ budget) : 0 ≤ D*(1+c) ∧ D*(1+c) ≤ budget := by
  constructor
  · exact mul_nonneg hD (by linarith)
  · have h := mul_le_mul_of_nonneg_left (show 1+c ≤ 2 by linarith) hD
    linarith

theorem gradientLoop (a b c d : ℚ) : (a-b)+(b-c)+(c-d)+(d-a) = 0 := by ring

#print axioms covarianceHalf
#print axioms lehmannPositive
#print axioms lehmannPair
#print axioms nonzeroModeNumerator
#print axioms brownianNormalization
#print axioms factorialStep
#print axioms thermalCoefficient
#print axioms signedFluxError
#print axioms gradientLoop

end AdiabaticPressure
