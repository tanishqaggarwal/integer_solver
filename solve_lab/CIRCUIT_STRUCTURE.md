# What the instance is, as an arithmetic circuit

Fully reverse-engineered and independently verified.

## Topology
- **383 combine gadgets** forming a binary tree over **384 leaf slots**
  (256 bit-selected + 128 hardwired zero pads). Root = gadget 149, whose two
  children hold 178 and 78 live leaves respectively.
- **256 boolean selector bits.** Each selector `b` gates exactly two pin atoms,
  forcing two free variables to hardcoded ~296-bit constants when `b = 1`
  (`b*(v - C) = 0`) and to zero when `b = 0` (`(1-b)*v = 0`).
- Each gadget is **fraction-free**, over register pairs (x1,y1),(x2,y2) -> (x3,y3):
  ```
  Rx = (x3+x1+x2+a2)*(x2-x1)^2 - (y2-y1)^2   == 0
  Ry = (y3+y1)*(x2-x1) - (y2-y1)*(x1-x3)     == 0
  ```
  gated by `sel = b1*b2` through 3 random combinations `sel*(c1*Rx + c2*Ry) = p*(free)`.
  Output is a 4-way MUX; the flags are OR-trees over the leaf bits.
- `a2 = x_24453` and `p = x_26064` are hardwired constant gates (see FINAL_CERTIFICATE.md).
- The root registers are compared against a hardcoded target pair; an atom forces
  OR(all selectors) = 1.

## The accumulator relation
The 256 leaf constants form a **self-combining chain** under the combine gadget:
`L_{i+1} = combine(L_i, L_i)`, verified for all 256 (0 mismatches). Consequently the
selected subset must accumulate, through the tree, to exactly the hardcoded target.
The selector bits are therefore the digits of a single 256-bit quantity, and the
system is satisfiable precisely when that quantity is the one the setter chose.

## Why every local method failed
The residual obstruction is the two-dimensional mismatch between the accumulated root
register pair and the target. It is invariant under every local move: no absorber,
Newton step, Grobner window, deflation, or selector-local edit can change it. This
retroactively explains the rank-1 "conserved obstruction", the reciprocal lock
r1*r2 = 1 mod p, the mod-p rigidity of all 766 slack absorbers (each exactly
`p*(free var)`), and why 39022 is the optimum of the local family.

## Soundness audit (all five checks)
1. **Division wires — weak form confirmed.** All 383 lambda-wires use
   `w*(x2-x1) = (y2-y1)`; there are ZERO wires of the strong form `w*(x2-x1) = 1`.
   When two operands coincide, both residuals vanish identically and the output
   registers are free. **This is the defect best_agentA_39022 already exploits.**
2. **Branch/flag manipulation — none.** All 256 selectors carry boolean atoms; all
   766 flag slots are gate-driven (638 with an exact 2-term equality atom, 128
   constant-0 pads); 0 free slots.
3. **Node degeneracy — unreachable.** Two children hold disjoint leaf subsets, so
   coincidence would need equal accumulations from disjoint sets. All 383 nodes fail
   even the magnitude test, and an exhaustive signed-digit search over all nodes and
   both orientations found 0 solutions.
4. **No curve-style consistency check exists** — off-curve intermediate registers pass
   every accumulator constraint. But backward solving from the target still terminates:
   all 1532 register slots are constrained (768 leaf wires + 764 internally linked),
   with **0 unconstrained slots**.
5. **Integer vs mod-p slack — sound.** All 1149 residual-gating constraints (383x3)
   carry a p-slack, 0 exceptions; 788 slacks are additionally forced to exact 0.

## Verified states
| file | score | note |
|---|---|---|
| best_agentA_39022.json | 39022/39033 | record; optimum of the local family |
| commonmode_39021.json | 39021/39033 | both gap atoms closed exactly over Z, core degenerate; residue = 2 atoms |
| nondegen_gapclosed_39013.json | 39013/39033 | gaps closed, non-degenerate core |

At best_agentA exactly **2 of 38133 atoms** are nonzero:
`7376877*x_642 + x_2099 - x_7068` and `x_4432 - x_19964 - x_28730`, both nonzero mod p.
