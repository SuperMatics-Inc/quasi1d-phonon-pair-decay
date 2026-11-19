import Mathlib.Data.Rat.Defs
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

namespace AdiabaticLinear3

theorem hermitianFactor (x y : ℚ) : ((x+y)/2)^2 - ((x-y)/2)^2 = x*y := by ring

theorem pressureBudget (t U : ℚ) : 2*t+2*U*2 = 2*(t+2*U) := by ring

theorem medianRange (p n z L : ℚ) (hL : L = p+n+z) (hp : 2*p ≤ L) (hn : 2*n ≤ L) :
    p-n ≤ z ∧ n-p ≤ z := by constructor <;> linarith

theorem medianSquare (d z : ℚ) (hd : 0 ≤ d) (hz : d ≤ z) : d^2 ≤ z*d := by
  have h := mul_le_mul_of_nonneg_right hz hd
  nlinarith

theorem chargeCost (a D cost : ℚ) (h : cost ≤ a*D/2) : -a*D+cost ≤ -a*D/2 := by
  linarith

theorem windingDamping (A m cost : ℚ) (hA : 0 ≤ A) (hm : 1 ≤ m) (h : cost ≤ A/4) :
    cost - m*A ≤ -m*A/2 := by
  have hmul : A ≤ m*A := by simpa using mul_le_mul_of_nonneg_right hm hA
  linarith

theorem zeroFluxMargin : (1 : ℚ) - 4*(1/9)/(1-1/9) = 1/2 := by norm_num

theorem momentCoefficients (m : ℚ) :
    m^2 - 3*(m-1)^2 + 3*(m-2)^2 - (m-3)^2 = 0 := by ring

theorem kineticAxisAndTwist : (1/4 : ℚ)/2/2 = 1/16 ∧ (2 : ℚ)^2 = 4 := by norm_num

#print axioms hermitianFactor
#print axioms pressureBudget
#print axioms medianRange
#print axioms medianSquare
#print axioms chargeCost
#print axioms windingDamping
#print axioms zeroFluxMargin
#print axioms momentCoefficients
#print axioms kineticAxisAndTwist
end AdiabaticLinear3
