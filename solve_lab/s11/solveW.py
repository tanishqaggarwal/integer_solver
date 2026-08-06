import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
P=L.P
p688=L.polys[688]; p406=L.polys[40608]
C688=p688[()]; k=p688[('x18956' and 18956,)] if False else None
# extract by var ids
def g(d,key): return d.get(key,0)
c0=g(p688,()); cX=g(p688,(14257,)); cY=g(p688,(18956,))
print("a688: c0=%d cX=%d cY=%d"%(c0,cX,cY))
q0=g(p406,()); qX=g(p406,(14257,)); qY=g(p406,(18956,))
qXX=g(p406,(14257,14257)); qXY=g(p406,(14257,18956)); qYY=g(p406,(18956,18956))
print("quad: XX=%d XY=%d YY=%d"%(qXX,qXY,qYY))
m=8863713
print("XY == -2m ?", qXY==-2*m, "  YY == m^2 ?", qYY==m*m)
print("qY == -m*qX ?", qY==-m*qX)
# W = X - m*Y  ; a688 = -W + c0  -> W = c0
W = c0
print("W required =", W)
# a40608 = W^2 + qX*W + q0 ?
val = W*W + qX*W + q0
print("a40608 at W:", val)
print("ZERO!" if val==0 else "NONZERO")
# roots of t^2+qX*t+q0
import math
D=qX*qX-4*q0
r=math.isqrt(D) if D>=0 else -1
print("disc is perfect square:", D>=0 and r*r==D)
if D>=0 and r*r==D:
    print("roots:", (-qX+r)//2, (-qX-r)//2)
