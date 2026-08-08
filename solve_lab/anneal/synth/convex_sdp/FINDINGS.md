# Convex (SDP / SOS / moment) relaxation of the encoder's gadgets

Shor/Lasserre-level-1 SDP `min <C,X> s.t. X>=0, diag(X)=1` solved in pure numpy
(mixing method; k>sqrt(2N) reaches the global SDP optimum). Looseness is rigorous
(mixing returns a feasible X, so <C,X> < E_min certifies a gap). Every gadget has
E_min=0 so gap = -SDP_opt. Variable-fixing gated on full ground-state enumeration.

## 1. Gadget-level exactness
| gadget | SDP_opt | gap | verdict |
|---|---|---|---|
| AND z=a*b | -0.125 | +0.125 | **LOOSE** |
| full-adder / 3:2 compressor (a+b+c-s-2d)^2 | 0 | 0 | **TIGHT** |
| not_equal (w=3,4) | 0 | 0 | **TIGHT** |

The AND is not degree-2 exact: bare diag=1 SDP omits the McCormick/RLT cuts
0<=z<=a, z<=b, z>=a+b-1 the Boolean-product hull needs (gap >= 1/8, explicit X*).
FA/compressor/not_equal are exact PROVABLY: in qubo.py every non-AND penalty comes
from add_square = (affine)^2 (verified for _ripple_eq, _wallace, _compress,
_wallace_eq, assert_zero, hence not_equal/mul_word/mul_eq/congruent). A sum of squares
of affine forms is a degree-2 SOS certificate => Lasserre-1 exact, gap 0 a priori. The
ONLY non-square penalty in the whole compiler is the AND (Rosenberg linearization).

## 2. modmul SDP bound vs exact (wallace)
| s | p | vars | #gs | SDP_opt | gap | removable(verified) |
|---|---|---|---|---|---|---|
| 4 | 11 | 73 | 18 | -24 | 24 | 1 |
| 6 | 41 | 153 | 95 | -54 | 54 | 1 |
| 8 | 163 | 255 | 403 | -96 | 96 | 1 |
Exact for no size: integrality gap never 1, GROWS with size (24->54->96, tracking AND
count). Optimal moment matrix deeply fractional (rank 15-25). The modmul is the needle
core; the degree-2 SDP does not see it.

## 3. SDP-persistent variables
1 removable at every realistic size (a single always-zero top compressor bit). At s=4
the SDP misses 9 of 10 truly-pinned bits (loose optimum washes out correlations). Sound
but weak detector. Reduction delivered: essentially none (1 qubit).

## 4. Convex-hull tightness
Arithmetic part: compact exact description (each add_square is one hyperplane; level-1
theta-body = true optimum). Single AND: hull is the integral McCormick simplex but needs
LINEAR RLT cuts the bare SDP omits. Product core: once thousands of ANDs share operand
bits the joint hull is the Boolean-quadric/correlation polytope (exponentially many
facets, not captured at fixed Lasserre level for growing size) -- the growing gap in §2
is the direct measurement.

## Verdict
Convex relaxations BOUND but do not SHRINK the problem. The arithmetic layer is convexly
exact (SOS) but was never the hard part; ALL hardness and the entire integrality gap
concentrate in the Boolean products the ANDs introduce -- exactly what the degree-2 SDP
cannot see. No subproblem of the F~8 frontier becomes convexly solvable in a way that
reduces annealed variables. Sound persistency removes 1 qubit.
