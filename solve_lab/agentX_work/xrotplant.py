import json,sys
d=json.load(open('xdata.json')); p=int(d['p']); A_=int(d['a']); N=int(d['N'])
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
S=[int(v) for v in sys.argv[1].split(',')]
P=None
for i in S: P=add(P,lad[i])
open(sys.argv[2],'w').write('%d %d\n'%P + ''.join('%d %d\n'%(x,y) for x,y in lad))
def f(S,j): return sum(1 for i in S if ((i-j)%256)<128)
print('planted |S|=%d  S=%s'%(len(S),S))
print('balanced rotations:',[j for j in range(128) if f(S,j)==5])
