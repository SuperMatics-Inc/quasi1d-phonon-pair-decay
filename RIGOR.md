# Exact charge-gauge bounds on finite-width cylinders

**PROVED.** For the positive-mass spinful Holstein and bond SSH models on a periodic axial cylinder, the canonical Yang pair density matrix has an explicit norm bound uniform in cylinder length and phonon energy at any fixed positive temperature and finite transverse width. The equilibrium axial helicity has a finite-length exponential upper bound in absolute value. Thus the exact positive-temperature pair-ordering threshold on a chain or any fixed-width ladder is zero, for all positive `Ω`, including at fixed interior density and fixed nonzero `t,U`. The flux result extends to every family whose aspect ratio tends to infinity. These statements give no isotropic two-dimensional improvement. The pressure cost in the complex-gauge argument does not become small with `Ω` in general, as shown by an exact Holstein lower bound.

## 1. Setup and statements

Let columns `j∈Z/LZ`, `L≥3`, each have `W≥1` sites. There is one axial bond per site joining the corresponding site of the next column (`WL` axial bonds). All other links are within a column. They can be the usual open or periodic transverse links. Take spin-degenerate hopping `−t` with `t>0` and independent neutral harmonic oscillators of positive energy `Ω`. Oscillator energy coordinates `b=FX` have coupling scale `U=F²/K>0`, where `K` is the spring constant. Holstein uses `b_i n_i` on sites; bond SSH uses `b_e B_e` on undirected edges. The electronic hopping bilinear on an unoriented link at real Peierls phase `η` is `B_e(η)=Σσ(e^(iη)c†_(jσ)c_(iσ)+e^(−iη)c†_(iσ)c_(jσ))`, with `‖B_e(η)‖≤2`; axial bonds are oriented consistently. There are no direct electron interactions. Set `τ=k_B T>0`, `β=1/τ`, and `E=t` (Holstein) or `E=t+2U` (SSH). Divide out free oscillator traces; additive constants in the free energy do not affect flux derivatives.

Fix any canonical particle number `N`. Let `Γ` be the electronic pair density matrix indexed by ordered spin-orbital pairs `i<j`, with entries `⟨c_i†c_j†c_l c_k⟩` for `k<l`. Let `F_N(θ)=-τ log Z_N(θ)` for uniform real electron phase `θ/L` per axial link, no transverse flux. Put

`a=min{1,1/(8βWE)}`, `δ=exp(−aL/2)`.

Then at zero flux, or the same real flux in the pair estimate,

`0≤Γ`,  `‖Γ‖ ≤ 8 W² [(1+exp(−a/2))/(1−exp(−a/2))]²`.   (1)

For `δ≤1/9`,

`|F_N''(0)| ≤ 4τ δ(1+δ)/(1−δ)^3`,

`|ρ_(s,x;L,W)| ≤ τ (L/W) δ(1+δ)/(1−δ)^3`.                (2)

In this normalization `ρ_(s,x;L,W)=(L/(4W))F_N''(0)`: the pair phase changes by twice the electron phase. Equation (1) excludes an extensive Yang eigenvalue at fixed `W` as `L→∞`, uniformly in `N` and `Ω>0`. It also excludes one at fixed `t,U,τ` along `W=o(L^(1/3))`: `a^−1=O(W)` and `‖Γ‖=O(W^4)` while volume is `LW`. Equation (2) forces axial helicity to zero whenever `L/W→∞`, since `a≥min{1,1/(8βE)}/W`. The same proof applies in a fixed-μ grand-canonical trace, but there is no ensemble conversion in the fixed-density claims. The convention of assigning `T_pair=0` means there is no positive-temperature pair order, not that the ground state has it.

## 2. Exact cost of the Hermitian part

At fixed real phases and for any selection of axial bonds, multiply the selected full hopping coefficients by independent real `s_e≥1`; for SSH this multiplies both `−t` and the force coefficient on that link. Let `Z(s)` denote the normalized partition function. For each scaled link,

`∂_(s_e) log Z(s) ≤ 2βt + 4β U s_e` (SSH), or `≤2βt` (Holstein).   (3)

The bare part follows from the norm bound. In SSH the free normalized harmonic Gaussian covariance for the energy coordinate is the periodic kernel

`D(u,v) = (U/β) Σ_(m∈Z) Ω²/((2πm/β)²+Ω²) exp(2π i m(u−v)/β)`.

It is nonnegative pointwise: the equally normalized kernel is `(UΩ/2)[exp(−Ωr)+exp(−Ω(β−r))]/[1−exp(−βΩ)]` for `r=|u−v|∈[0,β]`; `∫_0^β D(u,v)dv=U`. This follows by solving the periodic harmonic Green function or summing its Fourier coefficients. Gaussian integration by parts in the harmonic trace gives, at fixed real parameters,

`∂_(s_e) log Z(s) = βt⟨B_e⟩ + s_e ∫_0^β ∫_0^β D(u,v) ⟨T B_e(u)B_e(v)⟩ du dv`.

To check the second factor directly, write the oscillator-dependent electronic evolution as its Dyson series with coupling `−s_e b_e(u) B_e` in the exponent. The derivative in `s_e` inserts `−b_e(u)B_e`; differentiating the evolution with respect to `b_e(v)` inserts `−s_e B_e`. The two negatives cancel and there is no `1/2`, because the time slots of the newly contracted pair are both integrated. Gaussian exponential moments of `∫|b|` at finite volume dominate the source derivatives. The conditional time-ordered functions average to the full Gibbs trace functions by inserting the bounded electronic operators in the same harmonic trace. The fixed-`N` sector is unchanged since every `B_e` preserves it. For `r∈[0,β]`, Schatten Hölder bounds `|Tr(e^(−(β−r)H)B_e e^(−rH)B_e)|/Z≤‖B_e‖²≤4`. These facts give (3), including after other bonds were rescaled. For Holstein the scaled axial bonds do not change any site oscillator force, so the second term is absent.

Set `s_e=cosh x_e`, with `|x_e|≤1`. Integrating (3) shows the one-bond cost in `log Z(s)−log Z(1)` is bounded by `2β t(cosh x_e−1)+2β U sinh² x_e`, with the `U` term omitted for Holstein. For `|x|≤1`, `cosh x−1≤x²` and `sinh² x≤2x²`: the series give the bounds with constants `cosh 1−1` and `sinh² 1`, respectively; using `2<e<3` gives `cosh1<7/4` and `sinh1<4/3`. Thus

`log Z(s)−log Z(1) ≤ 2β E Σ_e x_e²`.                     (4)

This uses neither a bound on oscillator coordinates nor a trace-norm estimate on the SSH hopping operator.

## 3. Complex gauge trace inequality and oscillator cutoff removal

For finite-dimensional `A=A†`, `C=C†`, use the Lie product formula and Schatten Hölder:

`|Tr O e^(−β(A+i C))| ≤ ‖O‖ ‖e^(−β(A+i C))‖_1 ≤ ‖O‖ Tr e^(−β A)`.   (5)

Indeed for `n` factors `e^(−β A/n)e^(−iβ C/n)` their Schatten `n` norms equal `(Tr e^(−β A))^(1/n)`; the unitary factor has no effect. Take the product limit. With physical electron phase `η` and imaginary part `x_e`, the two oriented Peierls coefficients are `e^(iη−x_e)` and `e^(−iη+x_e)`. Their Hermitian part on the bond is exactly `cosh x_e` times the real-phase bond. Multiplying by a self-adjoint oscillator coordinate does not change the algebra because it commutes with the electronic operators.

Apply (5) first to oscillator Fock truncations on all sites or links, retaining the first `R+1` levels of each oscillator. No complex-flux limit is needed. For real phases and bounded `s_e`, full and compressed Hamiltonians are bounded below as forms by `A_0/2−C` with `A_0=Ω Σ a†a`, `C<∞` at fixed volume and fixed parameters. To see this use `‖a ψ‖²=⟨ψ,a†aψ⟩` and the bounded electronic coefficients to estimate linear oscillator forms by `ε⟨A_0⟩+C_ε‖ψ‖²`. Thus the full real Hamiltonian has compact resolvent and a finite heat trace, and the finite Fock subspaces increase to a form core. By the min-max principle the ordered compressed eigenvalues decrease to the full eigenvalues and are bounded below by the ordered eigenvalues of `A_0/2−C`. Convergence of the positive heat traces follows by summable domination, and any compressed trace with identity is at most its full real counterpart (before oscillator normalization). Dividing both sides at each step by the **same full free oscillator trace** preserves the inequalities. For bounded electronic insertions also use trace-norm convergence: the Ritz resolvent is the orthogonal projection of the full resolvent solution in the shifted form norm onto the finite subspace, hence converges strongly by the form core. The padded heat semigroup then converges strongly by continuous functional calculus of the bounded resolvents, with the function value at zero set to zero. Positive strong convergence plus convergence of traces gives trace-norm convergence: bound each off-diagonal finite-rank/complement block in trace norm by the geometric mean of the two diagonal traces via Schatten Cauchy–Schwarz, and take a finite projection approximating the limiting trace. This covers the insertion and the adjoint even if they do not commute with one another. All constants in (4)–(5) are independent of the Fock cutoff after this limit. The untruncated real partition function is continuous in real phases, as also follows by the finite-volume Gaussian/Dyson representation dominated by the oscillator exponential moments.

## 4. Pair matrix: axial charge transport

Consider the monomial `O=c_i† c_j† c_l c_k`, of norm at most one, and assign an integer `q_x` to each column: the number of creators minus annihilators in that column. `Σ q_x=0`. In a finite truncation use the invertible similarity `S=exp(Σ_x φ_x n_x)` on the electrons, where `n_x` is the spin-summed occupation of the column and `φ_x` is real. Then `S O S^−1=e^(Σ_x q_x φ_x) O`. Along an axial edge the forward scalar changes by `e^(φ_(x+1)−φ_x)`, the backward scalar by its reciprocal. Horizontal differences determine the Hermitian part in (5), while within-column links and oscillator-only terms are unchanged. Cyclic trace, (5), real cutoff removal, then (4) imply

`|⟨O⟩| ≤ exp[Σ_x q_x φ_x + 2β W E Σ_x (φ_x−φ_(x+1))²]`   (6)

when all horizontal differences have magnitude at most one. The signs are fixed by applying the similarity to both the insertion and the Hamiltonian, not just the Hamiltonian.

Let `q_x=Q_x−Q_(x−1)` with integer cyclic `Q_x`; for instance form partial sums of `q`. Choose an integer median `m` of the `Q_x`, and set `D=Σ_x |Q_x−m|`. Let `p,n,z` count entries strictly above, strictly below, and equal to `m`. On the first two sets put `g_x=+1,−1`; on the median set put `g_x=−(p−n)/z` if `z>0`. Then `Σ g_x=0`, `|g_x|≤1` because `|p−n|≤z` follows from `p,n≤L/2`, and `Σ Q_x g_x=D`. If `z=0` the two half sets have equal count. Also `Σ g_x²≤p+n+|p−n|≤2D`, the last step using integrality. Choose `φ_x−φ_(x+1)=−a g_x`, which is periodic because the differences sum to zero. Summation by parts gives `Σ q_x φ_x=−aD`. The pressure cost is at most `4β W E a²D≤aD/2`. Equation (6) gives `|⟨O⟩|≤e^(−aD/2)`.

Here `D` is precisely the minimum, over the two bijections between the creation and annihilation labels, of the sum of the circular column distances. Any matching creates a net integer edge flow of divergence `q`, with total absolute flow at most its path length. Every edge flow with these divergences has form `Q_x+b` for an additive constant, minimized in `ℓ¹` by `b=−m`. Conversely, route `Q_x−m` in its sign direction. Its nonzero directed edges cannot form a cycle, because a directed cycle on the ring requires all edges to be strictly of one sign, contradicting the median. Decompose this integer flow into unit paths from excess to deficit by starting at an excess, following occupied directions, and removing the path when it reaches a deficit. At each step flow conservation at nondeficits ensures continuation; no directed cycle ensures termination. First match opposite labels in the same column at zero cost. The paths now give a matching with sum of circular distances at most their length `D`. This proves the equality, including coincident labels and the case `q=0`.

For each creator row, sum the bound over both ordered annihilator labels and bound `exp(−a min(d_1,d_2)/2)≤exp(−a d_1/2)+exp(−a d_2/2)` using the two total matching distances. Each column has `2W` spin orbitals and the circle distance counts at most two columns at each positive distance. Thus the row sum is at most `2(2W)^2[(1+e^(−a/2))/(1−e^(−a/2))]^2`. Hermiticity and the symmetric row-sum bound give the operator norm bound (1); positivity is the usual Gram property of the two-electron annihilation operators. This excludes delocalized pair wavefunctions as well as fixed relative-coordinate ones. By the same construction, fixed-support fields of opposite nonzero charge in separated column intervals have exponentially decaying correlations: a ramp across the separating gap with constant field on the supports and compensating return ramp on the opposite arc produces the analogue of (6) with a cost quadratic in the ramp's slope. Alternatively expand each bounded field in finitely many charge monomials in their supports; the matching cost is at least the absolute total charge times the nearest cyclic column separation, regardless of neutral pairs within each support. This supplies exponential decay, not just absence of an extensive Yang eigenvalue at fixed width.

## 5. Exponential suppression of axial flux response

At finite Fock cutoff the real `θ` partition function is `2π` periodic, because changing `θ` by `2π` is a shift by the unitary axial gauge of phases `2π j/L`. The complex continuation is entire as a finite-matrix exponential, and the same similarity makes it `2π` periodic for complex `θ`. Let `z_m` denote the untruncated real partition function's Fourier coefficients in the convention `Z(θ)=Σ_m z_m e^(imθ)`. For finite cutoff use contour shift `θ→θ−i sign(m) aL` to bound the `m`th coefficient. On that contour (5) compares to the **real** Hamiltonian at the same real flux with all axial links scaled by `cosh a`. Its trace after cutoff is at most the untruncated scaled partition function. Divide by the full free oscillator trace and use (4) to get

`|z_m^R|≤Z_max exp(2β E WL a²−|m| aL)`, `Z_max=max_real θ Z(θ)`.

At finite real flux the compressed traces converge, bounded by `Z_max`, so dominated convergence yields the same bound on `z_m`. For `|m|≥1` use `2β W E a≤1/4` to obtain the simpler `|z_m|≤Z_max δ^|m|`. Hence Fourier convergence is absolute, reproduces the continuous real partition function, and justifies termwise differentiation twice. In fact the flux is analytic on a nonzero strip by this same bound; no untruncated trace at complex flux is required.

The nonconstant flux deviation is `sup|Z(θ)−z_0|≤2Z_max δ/(1−δ)`. In particular `Z(0)≥Z_max[1−4δ/(1−δ)]≥Z_max/2` if `δ≤1/9`. Using the maximum at real flux in this step rather than a seam estimate avoids an additional factor exponential in `W`. The model is time reversal invariant at zero flux, so `Z'(0)=0` (complex conjugation reverses the real twist in the oscillator and electronic occupation bases). Finally

`|Z''(0)|≤2Z_max Σ_(m≥1) m² δ^m = 2Z_max δ(1+δ)/(1−δ)^3`.

The series identity follows by twice applying `δ d/dδ` to the geometric series. With `F_N''(0)=−τ Z''(0)/Z(0)`, this gives (2). No nonnegativity of the finite-volume curvature is needed.

## 6. What the square-lattice route would still need

**PROVED obstruction to the same pressure comparison on squares.** For Holstein on even squares, at `N=V/2`, `0<U≤t/4`, `0<τ≤t/64`, LINEAR2 gives `t⟨K_el⟩/V>t/4`, uniformly in `Ω`. Rotation separates half into the axial kinetic operator `K_x`. At zero real flux, call `Z_s` the real partition function with the axial hopping scaled by `s`. Its exact derivative at `s=1` is `∂_s log Z_s|_1=βt⟨K_x⟩`. Therefore

`lim_(α→0) [log Z_(cosh α)−log Z_1]/(V α²) = βt⟨K_x⟩/(2V) > βt/16`.

Thus an improved Hermitian-part **cost** vanishing with `Ω` cannot hold uniformly near `α=0` even at this fixed regular bare density. The quadratic cost is a comparison to real increased hopping, not the true flux curvature. The Schatten inequality discards the anti-Hermitian current information which can cancel the order-`t` diamagnetism. On `W=L` the proven Fourier and pair estimates involve no decreasing function of `Ω`; helicity need not vanish by this proof. **OPEN:** a square-lattice upper bound with exponent strictly above `1/2` at fixed regular density. The cylinder theorem is an exact low-dimensional regime, not an interpolation to a square and not a priority claim about the complex-gauge method.

## 7. A signed pressure comparison on arbitrary graphs, including squares

**PROVED comparison lemma.** Let the electronic Hilbert space be any finite canonical sector, with a Hermitian electronic Hamiltonian `H_e(A)` at any fixed real Peierls field `A`. For each of `M` independent oscillators let `U_j>0`, `Ω_j>0` and let its energy-coordinate `b_j` multiply a bounded Hermitian electronic operator `B_j(A)`. An arbitrary bounded Hermitian electronic source independent of the oscillator coordinates can be included in `H_e(A)`. For any chosen subset `S` of oscillators, write `F_qu` and `F_mix` for the free energies with these oscillators quantum or replaced by static Gaussians, respectively; the remaining modes are unchanged and quantum. Normalize every quantum trace by the corresponding free harmonic trace and every static coordinate by its free Gaussian integral. At every `τ=1/β>0`,

`0 ≤ F_mix−F_qu ≤ τ Σ_(j∈S) log[sinh(βΩ_j/2)/(βΩ_j/2)] ≤ Σ_(j∈S) min{Ω_j/2, Ω_j²/(24τ)}`.  (9)

This holds before the volume limit, uniformly in the electron number and real field. For identical energies on a square lattice and all modes replaced, the energy-density error is at most `min{Ω/2,Ω²/(24τ)}` per site for one Holstein oscillator per site and twice this for two SSH bonds per site. The proof has two distinct parts. The constant `1/24` does not come from bounding a particular electronic Duhamel function.

The bounds are in the lineage of Berezin's convex trace inequalities (Math. USSR-Sb. **17**, 269–277, 1972) and Lieb's comparison of classical and quantum free energies ([10.1007/BF01646493](https://doi.org/10.1007/BF01646493)). The upper trace bound below is the Golden–Thompson comparison with a normalized free coordinate integral ([10.1103/PhysRev.137.B1127](https://doi.org/10.1103/PhysRev.137.B1127); [10.1063/1.1704727](https://doi.org/10.1063/1.1704727)). The lower bound is also derived by a covariance interpolation to exhibit the sign of each nonzero frequency. No priority is claimed for the comparison.

### 7.1 Lower bound: change the nonzero Matsubara modes using a physical auxiliary model

For mode `j` put `D_0(u-v)=U_j/β`. Let `D_q` be the positive periodic quantum covariance in section 2, of zero-mode integral `U_j`. For `s∈[0,1]`, assign independently to each replaced mode a static variable of covariance `(1−s)D_0` and a harmonic quantum coordinate with covariance `sD_q`, and couple their sum to `B_j`. Equivalently use the original independent variables with prefactors `sqrt(1−s)` and `sqrt(s)`. At each fixed static vector this is a self-adjoint, time-independent augmented Hamiltonian, including all the quantum modes. In particular, the expectation in the Gaussian differentiation formula below is not taken in a Gibbs state for a single time-dependent phonon history.

For finite volume write `p(s)=log Z_s` for the normalized canonical partition. The path representation and independent Gaussian integration by parts, as in section 2, give on `0<s<1`

`p'(s) = (1/2) Σ_(j∈S) ∫_0^β du ∫_0^β dv (D_q−D_0)(u−v) G_(j,s)(u−v)`,

where `G_(j,s)` is the time-ordered two-point function of `B_j` in the full augmented state. The factor `1/2` follows from differentiating the two square-root prefactors; the signs from the two insertions cancel. Equivalently this is differentiation in covariance. Bounded `B_j` and Gaussian continuity and exponential moments justify the insertions, derivatives in the open interval, and continuity at the endpoints. Conditioned on static coordinates, `G` is a trace over a discrete full quantum spectrum. Its ensemble average uses the **positive conditional partition weights**. The static variables introduce no extra time dependence.

For `ν_m=2πm/β`, take the Fourier transform on one period with no prefactor on the integral. The covariance difference has coefficients `0` at `m=0` and `U_j Ω_j²/(ν_m²+Ω_j²)` otherwise, by the oscillator Green function. For `m≠0`, a conditional Lehmann expansion gives

`G_hat_(j,s)(ν_m) = Z^(-1) Σ_(r,l) |(B_j)_(rl)|² (e^(−βE_r)−e^(−βE_l))(E_l−E_r)/[(E_l−E_r)²+ν_m²] ≥ 0`.

Terms with equal energies are zero. The imaginary pieces cancel on pairing `(r,l)` and `(l,r)`, so the sum is over ordered pairs with no additional doubling. For example the bound on the absolute integrated ordered formula is `β ||B_j||²` by trace Hölder and permits summing before integration. Averaging over the fixed static coordinates preserves positivity. The decaying (`O(m^(-2))`) covariance coefficients allow integrating the Fourier series termwise, giving

`p'(s) = (1/2) Σ_(j∈S) Σ_(m≠0) [U_j Ω_j²/(ν_m²+Ω_j²)] G_hat_(j,s)(ν_m) ≥ 0`.

Thus `F_qu≤F_mix`. Comparing subsets by freezing one more oscillator or varying the two weights separately gives the same sign for each change. This also shows why a signed Duhamel kernel is useful for pressure, even though its sign does not carry over to a physical current-inserted second derivative in flux.

The same construction with two independent quantum copies of each varied oscillator proves monotonicity of the normalized pressure in each `Ω_j` at fixed `U_j`: the zero-mode coefficient is unchanged, while `Ω²/(ν_m²+Ω²)` increases strictly with positive `Ω` for each `m≠0`. Normalization matters for this monotonicity. If every varied `B_j` commutes with the other electronic operators and with the electronic Hamiltonian, the difference is identically zero: resolve their joint eigenbasis and complete the square in each oscillator or static integral.

### 7.2 Upper bound: heat trace with a flat Brownian bridge

Keep the unselected oscillators and the electrons in one Hilbert space. For unshifted Hamiltonians the selected coordinates have kinetic operator `K=−Σ_(j∈S)(U_j Ω_j²/2) ∂_(b_j)^2`, and at fixed `b` the potential operator is

`V(b)=Σ_(j∈S)b_j²/(2U_j)+H_other+Σ_(j∈S)b_j B_j`.

Here `H_other` contains `H_e` and all unselected oscillators and their couplings. First take a finite Fock cutoff for the unselected modes, giving a finite-dimensional matrix potential in the selected coordinates. The matrix heat kernel Feynman–Kac representation for `K+V` uses Brownian bridges with increments of diffusion variances `U_j Ω_j²`; the free diagonal kernel is `c=Π_(j∈S)(2πβU_jΩ_j²)^(-1/2)`. For a continuous path put `b(u)=b+ξ(u)` with `ξ` a zero-endpoint bridge. Trace Hölder on the time-ordered product and the arithmetic–geometric mean bound give a pointwise trace estimate

`|Tr T exp(−∫_0^β V(b+ξ(u))du)| ≤ (1/β)∫_0^β du Tr exp(−β V(b+ξ(u)))`.

Indeed, for `n` equal time slices the Schatten `n`-norm of the `k`-th matrix exponential is `(Tr exp(−βV(b_k)))^(1/n)`. The time-sliced ordered product converges along each continuous path. The heat-kernel formula follows by this ordering (or the semigroup Feynman–Kac formula for continuous, Hermitian matrix potentials bounded below); the quadratic confining potential dominates its linear bounded matrix terms. After integration over `b` the right-hand side is finite and independent of the bridge by translation invariance, so the trace bound is integrable and gives `Tr exp(−β(K+V)) ≤ c ∫ db Tr exp(−βV(b))`. The lower bound on all real full Hamiltonians comes from absorbing the bounded matrix coupling in the oscillator quadratic forms. The Fock form-core and monotone Ritz trace argument in section 3 removes the cutoff; for the mixed side the same lower bound leaves a Gaussian-integrable upper envelope in the static coordinates. Alternatively, the same unselected free factors could have been canceled at finite cutoff.

For each selected oscillator the free harmonic trace is `(2 sinh(βΩ_j/2))^(-1)` and its static Gaussian normalization integral is `sqrt(2π U_j/β)`. Thus the inequality for normalized partitions is `Z_qu≤Z_mix Π_(j∈S) 2sinh(βΩ_j/2)/(βΩ_j)`. In particular the upper bound does not assume a sign of any current correlation or of individual electronic histories. For `x≥0`, `sinh x/x=∫_0^1 cosh(xs)ds≤e^x`. Also `(2j)(2j+1)≥6j` for `j≥1` gives `(2n+1)!≥6^n n!`, hence `sinh x/x≤e^(x²/6)` term by term. Substitute `x=βΩ_j/2` and multiply by `τ` to get (9); the formula is continuous at `Ω_j=0`.

### 7.3 Why this does not extend the square-lattice transition bound

**PROVED information obstruction.** Even with the sign and uniformity in real flux in (9), an intensive pressure comparison cannot bound the `O(1)` square-torus curvature from a static reference with zero helicity. For any error envelope `ε>0` per site and constant `D>0`, consider real, even, periodic test free energies `F_mix(θ)=0` and `F_qu(θ)=−D(1+cos θ)`. For `V≥2D/ε` they obey `0≤F_mix−F_qu≤Vε` uniformly in `θ` and `F_qu''(0)=D`, independently of `ε`. These functions are a counterexample to the inference from pressure control alone, not model realizations or a counterexample to the physical conjecture. Gauge invariance itself does not repair this: the diagonal unitary `exp(i Σ_x φ_x n_x)` changes edge phases by `φ_x−φ_y`, leaving normalized partition functions unchanged even for the SSH phonon-modulated bonds. Its exact Hessian therefore annihilates gradients. The harmonic flux on the torus is not a gradient, as its circulation round the axial loop is `θ`; this Ward identity does not constrain the missing curvature. A proof applying (9) to BKT would need extra control of the **flux-dependent** interpolation difference, at a scale smaller than the `O(1)` twist cost in two dimensions. Positivity in 7.1 is for a two-insertion auto-correlation and does not give a sign or an upper bound on that difference. **OPEN:** a better-than-square-root exact adiabatic exponent for the transition on squares at fixed nonzero bare parameters.
