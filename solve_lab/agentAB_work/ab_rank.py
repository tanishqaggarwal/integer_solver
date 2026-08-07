#!/usr/bin/env python3
import math
from math import comb, log2
def Vol(W): return sum(comb(256,j) for j in range(W+1))
def region(B): return sum(comb(256,j) for j in range(B+1,257))
print("Break-even: largest B whose ball-covering proof still costs < 2^126.5 (= solve outright)")
prev=None
for B in range(255,100,-1):
    R=region(B)
    best=min(math.log2(max(1.0,R/Vol(W)))+log2(comb(256,W//2)) for W in range(2,200,2))
    if best>=126.5:
        print("  first B whose proof costs >= 2^126.5 :", B, " (cost 2^%.1f)"%best)
        print("  so the cheapest-to-prove nonvacuous ceiling is  w <= %d"%prev)
        break
    prev=B
print()
print("Ranked ledger: cost of the complement mechanism vs what it proves")
print("  W    cost         proves w <=   null mass of excluded region   equiv. sigma above mean 128")
for W in (6,8,10,12,14,16,20,24,30):
    c=log2(comb(256,W//2)); B=255-W
    q=region(B)/2**256
    print("  %2d   2^%5.1f      %3d           2^%7.1f                    %+.1f sigma"
          %(W,c,B,log2(q),(B+1-128)/8.0))
print()
print("For comparison, the SAME budget spent on the low-weight side:")
for W in (6,8,10,12,14,16,20,24,30):
    print("  W=%2d  cost 2^%5.1f  proves w >= %d   (null mass 2^%.1f)"
          %(W,log2(comb(256,W//2)),W+1,log2(sum(comb(256,j) for j in range(W+1))/2**256)))
