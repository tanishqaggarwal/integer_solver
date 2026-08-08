#!/usr/bin/env python3
"""embed.py -- physical-qubit estimates.

Logical variable count is not the binding hardware constraint; the largest clique
in the Hamiltonian is.  Every penalty here is a square of a linear form, and a
square over c terms makes a K_c.  Pegasus and Zephyr both embed a clique of side
~12m out of ~24m^2 qubits, so a K_c costs about c^2/4 physical qubits, i.e. each
logical variable in a K_c carries a chain of length ~c/4.
"""
import random, json, math
from subsetsum import build as ss_build
from resources import marginal_window

# Pegasus P_m: ~24m(m-1) qubits embed a clique of side ~12m-10, so a K_c costs
# about c/6 physical qubits per logical variable.  Calibrated on P16 (Advantage):
# K_182 fills all 5,760 qubits => 31.6 physical per logical = 182/5.8.
phys = lambda v, c: v * max(1.0, c / 6.0)
HW = 4400

print(f"{'encoding':>36} {'logical':>11} {'clique':>7} {'~physical':>13} {'vs 4400':>9}")
print("-" * 82)

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
rnd = random.Random(7)
logs = [rnd.randrange(N) for _ in range(256)]
tgt = rnd.randrange(N)
print("IF THE DISCRETE LOGS WERE KNOWN  (modular subset-sum over 256 relations)")
rows = []
for mode, ch, lbl in (('binary', 1024, 'no compression (fewest qubits)'),
                      ('binary', 64, 'chunk-compressed'),
                      ('wallace', 16, 'full-adder compressed')):
    Q, _ = ss_build(logs, tgt, N, mode=mode, chunk=ch)
    st = Q.stats()
    ph = phys(st['vars'], st['max_clique'])
    rows.append(ph)
    print(f"{lbl:>36} {st['vars']:11,d} {st['max_clique']:7d} {ph:13,.0f} {ph/HW:8,.0f}x")
print(f"{'--> invariant across encodings:':>36} {'':11} {'':7} {min(rows):13,.0f} {min(rows)/HW:8,.0f}x")

print()
print("AS THE INSTANCE ACTUALLY STANDS  (logs unknown -- the full comb ladder)")
for mode, w in (('binary', 9), ('wallace', 8)):
    v, c, jr, mc = marginal_window(256, w, mode, neq=True, want_clique=True)
    M = math.ceil(256 / w)
    ph = phys(M * v, mc)
    print(f"{'comb ladder, ' + mode + f', w={w}':>36} {M*v:11,d} {mc:7d} {ph:13,.0f} {ph/HW:8,.0f}x")

print("""
Both rows of the first block land at ~1e5 physical qubits however you trade
cliques against ancillas -- shrinking the clique inflates the ancillas by the
same factor.  That is the real floor of the "dream" encoding.  The second block
is the same decision problem without the logs, and it is ~100x worse again.""")
