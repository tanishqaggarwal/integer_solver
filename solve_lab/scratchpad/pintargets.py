import sys, os
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import heal_harness as H
import atomlib as A
p = H.p

# the 5 ancestors of x_2099, x_19964
for x in [2099, 19964]:
    print(f"x_{x} ancestors: {sorted(H.anc.get(x,set()))}")

# What are the boolean/bit vars? load all_bits.json
import json
bits = set(json.load(open('all_bits.json'))) if os.path.exists('all_bits.json') else set()
print(f"total bit vars: {len(bits)}")

# For each pin target, how many ancestors are bits vs non-bits
for name,x in [('x_1308',1308),('x_23927',23927),('x_17601',17601),('x_19083',19083)]:
    a = H.anc.get(x,set())
    ab = a & bits
    anb = a - bits
    print(f"\n{name}: {len(a)} anc = {len(ab)} bits + {len(anb)} non-bits")
    print(f"   non-bit ancestors: {sorted(anb)[:40]}")

# Q1 = x_1308 - x_23927, Q2 = x_17601 - x_19083. Which free inputs do they jointly depend on?
a1,a2 = H.anc.get(1308,set()), H.anc.get(23927,set())
b1,b2 = H.anc.get(17601,set()), H.anc.get(19083,set())
print(f"\nQ1 (1308 vs 23927): union anc = {len(a1|a2)}, shared={len(a1&a2)}")
print(f"Q2 (17601 vs 19083): union anc = {len(b1|b2)}, shared={len(b1&b2)}")
# do Q1 and Q2 share free inputs?
print(f"Q1 ∩ Q2 free inputs: {sorted((a1|a2)&(b1|b2))[:30]}")
print(f"|Q1∪Q2 all free inputs| = {len((a1|a2|b1|b2))}")
