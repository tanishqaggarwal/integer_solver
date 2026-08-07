import json,os,sys
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
S=[int(v) for v in sys.argv[1].split(',')]
P=None
for i in S: P=add(P,lad[i])
out=sys.argv[2]
with open(out,'w') as f:
    f.write('%d %d\n'%(P[0],P[1]))
    for x,y in lad: f.write('%d %d\n'%(x,y))
print('planted S=%s |S|=%d -> %s'%(S,len(S),out))
