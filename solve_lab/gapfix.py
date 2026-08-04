#!/usr/bin/env python3
import heal_harness as H
import pickle
p=H.p
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.val[14853]=vA[14853]; H.val[31339]=vA[31339]
H.forward(); v=H.val
# gap 24105 = 4261533*x_6858 + x_25295 - 4261533*x_31339 ; x_25295 = x_31999*x_11559
print(f"x_31999 = {v[31999]}   (==p? {v[31999]==p}, wire? {31999 in wire}, %p={v[31999]%p if v[31999]!=p else 0})")
print(f"x_11360 = {v[11360]}   (==p? {v[11360]==p}, wire? {11360 in wire})")
tgt24105 = 4261533*(v[31339]-v[6858])   # need x_25295 = this
print(f"\ngap 24105: need x_25295 = 4261533*(x_31339-x_6858) = {tgt24105}")
print(f"  = x_31999 * x_11559 ; x_11559 = tgt/x_31999 ; divisible by x_31999? {tgt24105 % v[31999]==0}")
tgt27902 = 12846437*(v[14853]-v[1308])
print(f"gap 27902: need x_29967 = 12846437*(x_14853-x_1308) = {tgt27902}")
print(f"  = x_11360 * x_30163 ; divisible by x_11360? {tgt27902 % v[11360]==0}")
# Try: set x_11559, x_30163 to fix gaps (if divisible), forward, check
F0=set(H.fails())
if tgt24105 % v[31999]==0 and tgt27902 % v[11360]==0:
    H.val[11559]=tgt24105//v[31999]
    H.val[30163]=tgt27902//v[11360]
    H.forward(); F1=set(H.fails())
    print(f"\nAFTER activating gap slacks: {len(F1)} fails (was {len(F0)}); fixed {sorted(F0-F1)}, broke {sorted(F1-F0)[:10]}")
else:
    print("\nNOT cleanly divisible -> gap slacks p-granular or need lift")
