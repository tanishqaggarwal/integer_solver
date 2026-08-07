import json,os,numpy as np,random,itertools
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'xdata.json')))
p=int(d['p']); A_=int(d['a'])
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
    x3=(l*l-x1-x2)%p; return (x3,(l*(x1-x3)-y1)%p)
M=(1<<64)-1
k=np.fromfile('tbl3.bin',dtype=np.uint64)
print('table entries',len(k))
# block 1
ok=all(int(k[m])==lad[m][0]&M for m in range(256)); print('size-1 block 256/256:',ok)
# block 2 : full
off=256; bad=0; idx=0
for i in range(256):
    for j in range(i+1,256):
        if int(k[off+idx])!=add(lad[i],lad[j])[0]&M: bad+=1
        idx+=1
print('size-2 block %d entries, mismatches %d'%(idx,bad))
# block 3 : random sample
off3=256+32640
trip=list(itertools.combinations(range(256),3))
print('size-3 count',len(trip),'file remaining',len(k)-off3)
random.seed(1); bad=0
for _ in range(400):
    t=random.randrange(len(trip)); i,j,l=trip[t]
    if int(k[off3+t])!=add(add(lad[i],lad[j]),lad[l])[0]&M: bad+=1
print('size-3 random 400 mismatches',bad)
