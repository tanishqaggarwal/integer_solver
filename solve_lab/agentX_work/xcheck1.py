import json,os
HERE=os.path.dirname(os.path.abspath(__file__))
d=json.load(open(os.path.join(HERE,'xdata.json')))
p=int(d['p']); A_=int(d['a']); B_=int(d['b'])
lad=[(int(a),int(b)) for a,b in d['ladder']]; T=(int(d['T'][0]),int(d['T'][1]))
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
def neg(P):return (P[0],(-P[1])%p)
def limbs(v):return [ (v>>(64*i))&(2**64-1) for i in range(4)]
print('L5+L7 x limbs :',limbs(add(lad[5],lad[7])[0]))
r=add(T,neg(lad[3]))
print('T-L3  x limbs :',limbs(r[0]))
print('T-L3  y limbs :',limbs(r[1]))
