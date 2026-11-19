"""Exact rational guards for LINEAR3, plus exhaustive small cyclic matchings."""

import argparse
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import subprocess
import sys

MUTANTS = ('green-half', 'pressure-ssh', 'median-entry', 'median-balance', 'transport-sign',
           'winding-sign', 'fourier-denominator', 'fourier-moment', 'pair-phase',
           'kinetic-axis')


def check(condition, guard):
    if not condition:
        raise SystemExit(f'REJECT: {guard}')


def pressures(mutant):
    u, omega, t = Q(6, 7), Q(2, 5), Q(1, 3)
    for boltzmann in (Q(1, 10), Q(3, 4), Q(1, 2)):
        # Each reflected exponential integrates to (1-q)/omega.
        divisor = 1 if mutant == 'green-half' else 2
        g_sq = u*omega/divisor
        integral = g_sq*2*(1-boltzmann)/(omega*(1-boltzmann))
        check(integral == u, 'normalized periodic oscillator covariance')
    z = Q(3, 2)
    c = (z+1/z)/2
    s = (z-1/z)/2
    check(c*c-1 == s*s, 'reciprocal Peierls factors give a cosh Hermitian part')
    coefficient = 1 if mutant == 'pressure-ssh' else 2
    # Cost is 2 beta t(cosh-1) + 2 beta U sinh²; use cosh-1≤a² and sinh²≤2a².
    exact_envelope = 2*t + 2*u*2
    check(exact_envelope == 2*(t+coefficient*u), 'SSH quadratic pressure envelope')
    check(((3+Q(1, 2))/2-1) < 1 and (Q(4, 3))**2 < 2, 'rational hyperbolic coefficients')
    axial_divisor = 1 if mutant == 'kinetic-axis' else 2
    check(Q(1, 4)/axial_divisor/2 == Q(1, 16), 'axis and cosh factor in nonvanishing cost')


def circle(a, b, length):
    direct = abs(a-b)
    return min(direct, length-direct)


def transport(mutant):
    for length in range(3, 9):
        for i, j, k, ell in product(range(length), repeat=4):
            charge = [(x == i)+(x == j)-(x == k)-(x == ell) for x in range(length)]
            flow = [sum(charge[:x+1]) for x in range(length)]
            ordered = sorted(flow)
            med = ordered[(length-1)//2]
            if mutant == 'median-entry':
                med = (med+ordered[length//2])//2
            p = sum(x > med for x in flow)
            n = sum(x < med for x in flow)
            z = length-p-n
            check(z >= 1, 'chosen median is an entry of the flow')
            balance = 1 if mutant == 'median-balance' else -1
            g = [(x > med)-(x < med) if x != med else Q(balance*(p-n), z) for x in flow]
            check(sum(g) == 0 and all(abs(x) <= 1 for x in g), 'median correction gives a periodic ramp')
            d = sum(abs(x-med) for x in flow)
            matches = (circle(i, k, length)+circle(j, ell, length),
                       circle(i, ell, length)+circle(j, k, length))
            check(d == min(matches), 'circulating median flow versus independent pair matchings')
            check(sum(x*x for x in g) <= 2*d, 'quadratic bound on axial ramp')
            alpha = Q(1, 8)
            sign = -1 if mutant == 'transport-sign' else 1
            phi = [Q(0)]+[sign*alpha*sum(g[:x]) for x in range(1, length)]
            pairing = sum(x*y for x, y in zip(charge, phi))
            check(pairing == -alpha*d, 'cyclic summation by parts for the decaying gauge')
            check(pairing+2*alpha**2*sum(x*x for x in g) <= -alpha*d/2,
                  'charge decay after paying the pressure cost')


def flux(mutant):
    for length, a in ((3,Q(1)), (64,Q(1, 7)), (104,Q(1, 32))):
        for m in range(1, 7):
            sign = 1 if mutant == 'winding-sign' else -1
            check(a*length/4+sign*m*a*length <= -m*a*length/2,
                  'contour damping of nonzero winding harmonics')
    mult = 2 if mutant == 'fourier-denominator' else 4
    delta = Q(1, 9)
    check(1-mult*delta/(1-delta) == Q(1, 2), 'compare zero flux to the real maximum')
    power = 1 if mutant == 'fourier-moment' else 2
    for m in range(15):
        # Coefficients of (1-x)^3 sum_{j≥1} j² x^j are x+x².
        coefficient = sum(c*(m-k)**power for k, c in enumerate((1, -3, 3, -1)) if m-k >= 1)
        check(coefficient == int(m in (1, 2)), 'second flux derivative geometric series')
    pair = 1 if mutant == 'pair-phase' else 2
    length, width = Q(12), Q(3)
    check(width*length*(pair/length)**2 == 4*width/length, 'pair-phase stiffness normalization')


def run(mutant):
    pressures(mutant)
    if not mutant or mutant in ('median-entry', 'median-balance', 'transport-sign'):
        transport(mutant)
    flux(mutant)
    print('PASS: exact pressure, cyclic charge transport, and axial winding coefficients')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--mutant', choices=MUTANTS)
    mode.add_argument('--all-mutants', action='store_true')
    args = parser.parse_args()
    run(args.mutant)
    if args.all_mutants:
        for name in MUTANTS:
            out = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()), '--mutant', name],
                                 text=True, capture_output=True)
            check(out.returncode == 1 and out.stderr.startswith('REJECT:'), f'{name} missed intended guard')
            print(f'{name} exit 1: {out.stderr.strip()}')
