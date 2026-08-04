# Core Reduction — MAJOR BREAKTHROUGH

## The entire 39,033-equation trapdoor reduces to THREE scalar conditions.

forward_construct.py satisfies all wiring (39,013/39,033). The 20 remaining are
verifier squares E²=0. **All 20 are exact integer combinations of just three monsters:**

- **M1** = x_15298·x_11150 + x_4007  (≈ x_11150, the load L1)
- **M2** = x_15298·x_25739 − 6672769·x_29804  (≈ x_25739, load L2)
- **M3** = 537773·(x_15298·x_37758) − x_35605  (≈ 537773·x_37758, load L3)

(x_15298=1 in quadrant (1,1); C1/C2 pairs x_24468−x_13682, x_18956−x_37892 cancel.)

Since eq 22093 = M2, eq 25539 = 7·M3, eq 2071 = M1+5·M2, the core forces **M1=M2=M3=0**.

## The three loads are wiring-defined by two base gates S=x_35389, T=x_6671:
- L1 = x_11150 = 8646263·S + 1073965·T
- L2 = x_25739 = 10159099·S + 6926539·T
- L3 = x_37758 = 8272701·S + 5921311·T

All three ≡0 mod p ⟺ (2×2 minors invertible) ⟺ **S≡0 mod p AND T≡0 mod p**.

## S,T residues collapse to TWO deep control gates:
- S = x_33469·x_29356 − x_3558²,  T = x_32680 − x_11602 = x_27713·x_29322 − x_3558·x_1326
- x_29356 = x_29322², x_27762 = x_3558²
- So **S≡T≡0 mod p ⟺ x_3558≡0 mod p AND x_29322≡0 mod p**

Both control gates are DIFFERENCES with a free-input side:
- **x_29322 = x_14853 − x_12186**  (both FREE — set x_14853≡x_12186 mod p)
- **x_3558 = x_24908 − x_16742**    (x_16742 FREE — set x_16742≡x_24908 mod p)

## Quotient handles (private, only in core eqs) finish M1,M2,M3:
- M1=0: x_30317 = −L1/p   (needs L1≡0 mod p ✓)
- M3=0: x_2936 = 537773·L3/p  (needs L3≡0 mod p ✓)
- M2=0: x_5146 = L2/(6672769·p)  (needs L2≡0 mod **6672769·p** — extra mod-6672769 cond)

## STATUS
- Minimal residue fix (x_14853,x_16742 by <p) achieves S≡T≡0 mod p, M1=M3=0 exactly. VERIFIED.
- REMAINING obstacle 1: fixing x_14853,x_16742 breaks ~23 non-core wiring eqs (coupled through
  S,T cone). Heal requires S,T-cone handles + re-pin residues (alternation).
- REMAINING obstacle 2: M2 mod-6672769 condition (L2/p % 6672769 = 645924); tunable via the
  QUOTIENT parts of x_29322 (d) and x_3558 (m): 10159099·s'+6926539·t' ≡ 0 mod 6672769,
  where s'=S/p, t'=T/p are functions of d,m. 2 knobs, 1 quadratic mod-6672769 condition → solvable.

## Quadrant analysis (this session, continued)
The MUX control gates: x_7715 = OR(x_8599,x_21839), x_34554 = OR(x_25956,x_7304),
x_15298 = x_7715·x_34554. The override {24601:1,2081:1} forces x_15298=1 (quadrant 1,1).

- **Quadrant (1,1)** (x_15298=1): forward_construct → 39,013. Core = M1=M2=M3=0 ⟺ S,T≡0 mod p.
  Sparse-witness null-space solve (30, then 149, 902, 6114 inputs) is INCONSISTENT even across
  the full coupling closure — S,T residues are LINEARLY PINNED by the constraining wiring.
- **Quadrant (0,0)** (x_15298=0, all activators off): greedy → 39,006 (27 fail). The core loses
  its x_15298·load terms and reduces to a SMALLER gadget set:
    G1 = x_24468 − C1 − 12354891·x_34243,   G2 = x_18956 − C2 − x_32237,
    x_2300=0, x_9274=0, and a g_complex over {x_12186,x_14853,x_22162}.
  BUT Dixon/accumulate still ripples globally (27→118 fail) — same multi-role coupling.

## Conclusion on the core
The 20 verifier squares E²=0 are present in BOTH quadrants (different E). The control variables
(x_24468,x_18956,x_12186,x_14853,x_16742,…) are MULTI-ROLE: each appears in ~40 equations, so
setting them to satisfy the core ripples through the wiring. The sparse wiring-satisfying
solution is unique (rank 30), pinning the residues. No local move (up to the full closure)
reaches the core. This is the irreducible trapdoor: solving it appears to need either the
setter's witness or a global nonlinear solver beyond greedy/Dixon/null-space perturbation.
Best verified partial remains 39,013/39,033.
