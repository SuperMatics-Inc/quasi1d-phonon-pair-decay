"""Exact rational replay for RIGOR.md §7. No third-party dependencies."""

import argparse
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import subprocess
import sys


MUTANTS = {
    "sqrt-half": "covariance derivative from the two independent amplitudes",
    "lehmann-sign": "ordered Lehmann pair sign and magnitude",
    "lehmann-double": "ordered Lehmann pair sign and magnitude",
    "zero-mode": "unchanged zero mode and ordered nonzero modes",
    "mass-normalization": "Brownian to static free prefactor squared",
    "zero-point": "normalized quantum-to-static heat trace prefactor",
    "thermal-coefficient": "coefficient of the first nonconstant power",
    "flux-obstruction": "pressures agree per volume but curvatures differ",
    "gradient-loop": "gradient circulation is zero, unlike torus flux",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def replay(mutant):
    def test(name, condition):
        require(condition, MUTANTS[name])

    for a, b in [(Q(3, 5), Q(4, 5)), (Q(5, 13), Q(12, 13))]:
        a_prime, b_prime = -1 / (2 * a), 1 / (2 * b)
        if mutant == "sqrt-half":
            a_prime *= 2
        test("sqrt-half", a*a + b*b == 1 and a*a_prime == -Q(1, 2)
             and b*b_prime == Q(1, 2))

    for delta in (Q(1, 8), Q(3, 2), Q(9)):
        for nu in (Q(1, 2), Q(4)):
            # Any ordered Gibbs weights w0>w1 have this positive sign;
            # rational weights avoid certifying a rounded exponential.
            w0, w1 = Q(1), Q(1) / (1 + delta)
            direct = (w0 - w1) * delta / (delta**2 + nu**2)
            paired_real = direct - (w1-w0)*delta/(delta**2 + nu**2)
            ordered = direct
            if mutant == "lehmann-sign":
                ordered = -ordered
            if mutant == "lehmann-double":
                ordered *= 2
            test("lehmann-sign", 0 < ordered == paired_real/2)
            for u in (Q(1, 9), Q(5), Q(50)):
                cov = lambda omega: u*omega**2 / (nu**2 + omega**2)
                omega1, omega2 = Q(1, 5), Q(7, 4)
                observed_zero = u if mutant == "zero-mode" else u-u
                exact_increase = (u * nu**2 * (omega2**2-omega1**2)
                                  / ((nu**2+omega1**2)*(nu**2+omega2**2)))
                test("zero-mode", observed_zero == 0 < cov(omega1) < cov(omega2)
                     and cov(omega2)-cov(omega1) == exact_increase)

    # Factor 2*pi cancels between the squared Brownian prefactor and integral.
    for beta, u, omega in [(Q(1, 4), Q(2, 3), Q(5, 2)), (Q(9), Q(12), Q(1, 5))]:
        diffusion = u*omega**2 if mutant != "mass-normalization" else u*omega
        raw_ratio_sq = (u/beta) / (beta*diffusion)
        test("mass-normalization", raw_ratio_sq == 1/(beta**2 * omega**2))
        # Exact coefficient series in y=beta*omega for 2*sinh(y/2)/y.
        double = 1 if mutant == "zero-point" else 2
        for n in range(31):
            coeff = Q(double, 2**(2*n+1)*factorial(2*n+1))
            test("zero-point", n != 0 or coeff == 1)
            displayed = Q(1, 12) if mutant == "thermal-coefficient" else Q(1, 24)
            test("thermal-coefficient", n != 1 or coeff == displayed)
            require(factorial(2*n+1) >= (6**n)*factorial(n), "factorial recursion")
            require(coeff <= Q(1, 24**n * factorial(n)), "coefficient exponential majorant")

    for energy_error, curvature in [(Q(1, 43), Q(2)), (Q(5, 2), Q(1, 2))]:
        volume = int(2*curvature/energy_error) + 1
        for cosine in (Q(-1), Q(-1, 8), Q(1)):
            pressure_diff = curvature*(1+cosine)
            if mutant == "flux-obstruction":
                pressure_diff *= -1
            test("flux-obstruction", 0 <= pressure_diff <= volume*energy_error
                 and curvature > 0)
    for phases in [(Q(1), Q(2), Q(3)), (Q(-1, 2), Q(8), Q(5, 6), Q(3))]:
        gradient = [x-y for x, y in zip(phases, (*phases[1:], phases[0]))]
        if mutant == "gradient-loop":
            gradient[-1] = -gradient[-1]
        theta = Q(8, 5)
        harmonic = [theta/len(phases)]*len(phases)
        test("gradient-loop", sum(gradient) == 0 and sum(harmonic) == theta != 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--mutant", choices=MUTANTS)
    group.add_argument("--all-mutants", action="store_true")
    args = parser.parse_args()
    try:
        replay(args.mutant)
    except ValueError as error:
        print(f"REJECT: {error}")
        return 1
    print("PASS: pressure signs, heat trace normalization, and square-flux obstruction")
    if args.all_mutants:
        for mutant, message in MUTANTS.items():
            result = subprocess.run(
                [sys.executable, "-B", str(Path(__file__).resolve()), "--mutant", mutant],
                capture_output=True, text=True, check=False,
            )
            print(f"{mutant} exit {result.returncode}: {result.stdout.strip()}")
            if result.returncode != 1 or result.stdout.strip() != f"REJECT: {message}":
                return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
