import json,random,itertools
import numpy as np
d=json.load(open('xdata.json')); p=int(d['p']); A_=int(d['a'])
lad=[(int(a),int(b)) for a,b in d['ladder']]
def inv(z):return pow(z,p-2,p)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x1,y1=P;x2,y2=Q
    if (x1-x2)%p==0:
        if (y1+y2)%p==0:return None
        l=(3*x1*x1+A_)%p*inv(2*y1%p)%p
    else: l=(y2-y1)%p*inv((x2-x1)%p)%p
    x3=(l*l-x1-x2)%p;return (x3,(l*(x1-x3)-y1)%p)
def sp(s):  # signed ladder point
    m=s>>1; P=lad[m]
    return P if not (s&1) else (P[0],(-P[1])%p)
M=(1<<64)-1
k=np.fromfile('stbl.bin',dtype=np.uint64); print('signed table entries',len(k))
# a=1 block: 256 entries, s0 = 0,2,...,510
ok=all(int(k[i])==lad[i][0]&M for i in range(256)); print('a=1 block 256/256:',ok)
# a=2 block: enumerate in the same order
idx=0; bad=0; off=256
pairs=[]
for s0 in range(0,512,2):
    for s in range(((s0>>1)+1)<<1, 512):
        pairs.append((s0,s))
print('a=2 count',len(pairs))
for (s0,s) in pairs:
    if int(k[off+idx])!=add(sp(s0),sp(s))[0]&M: bad+=1
    idx+=1
print('a=2 block mismatches',bad)
# a=3 random sample
off3=256+len(pairs)
trip=[]
for s0 in range(0,512,2):
    for s1 in range(((s0>>1)+1)<<1,512):
        for s2 in range(((s1>>1)+1)<<1,512): trip.append((s0,s1,s2))
print('a=3 count',len(trip),' file remaining',len(k)-off3)
random.seed(7); bad=0
for _ in range(400):
    t=random.randrange(len(trip)); s0,s1,s2=trip[t]
    if int(k[off3+t])!=add(add(sp(s0),sp(s1)),sp(s2))[0]&M: bad+=1
print('a=3 random 400 mismatches',bad)
k.sort(); k.tofile('stbls.bin'); print('sorted -> stbls.bin')
