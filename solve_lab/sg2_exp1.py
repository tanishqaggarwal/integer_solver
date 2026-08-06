#!/usr/bin/env python3
"""Experiment 1: in the forward model, set x_7068:=x_2099 and x_4432:=x_19964,
partners x_17325=x_9413=0, forward(), check what fails."""
import heal_harness as H
p = H.p

vA = H.loadd('best_agentA_39022.json')
for v in H.freeinp:
    H.val[v] = vA.get(v, 0)
H.forward()
F0 = H.fails()
print(f"baseline: {len(F0)} fail: {F0}")

# read computed gate values
x2099 = H.val[2099]; x19964 = H.val[19964]
x7068 = H.val[7068]; x4432 = H.val[4432]
print(f"x_7068={x7068%p} (mod p)  x_2099={x2099%p} (mod p)  equal? {(x7068-x2099)%p==0}")
print(f"x_4432={x4432%p} (mod p)  x_19964={x19964%p} (mod p)  equal? {(x4432-x19964)%p==0}")
print(f"gap G1=x_7068-x_2099 = {x7068-x2099}")
print(f"   /p = {(x7068-x2099)//p}, %p = {(x7068-x2099)%p}, %7376877={(x7068-x2099)%7376877}")
print(f"gap G2=x_4432-x_19964 = {x4432-x19964}")
print(f"   /p = {(x4432-x19964)//p}, %p = {(x4432-x19964)%p}")

# Now set x_7068 := x_2099, x_4432 := x_19964, partners 0
print("\n--- setting x_7068:=x_2099, x_4432:=x_19964, partners=0 ---")
H.val[7068] = x2099
H.val[4432] = x19964
H.val[17325] = 0
H.val[9413] = 0
H.forward()
# did x_2099, x_19964 change?
print(f"after forward: x_2099={H.val[2099]} (changed? {H.val[2099]!=x2099})")
print(f"after forward: x_19964={H.val[19964]} (changed? {H.val[19964]!=x19964})")
print(f"after forward: x_7068={H.val[7068]}  x_4432={H.val[4432]}")
F1 = H.fails()
print(f"after: {len(F1)} fail: {F1[:40]}")
