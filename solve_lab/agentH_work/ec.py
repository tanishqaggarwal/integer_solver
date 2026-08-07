"""EC arithmetic on y^2 = x^3 + B over p, plus group-order determination."""
import json, math
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
K=97553848499418123410591666447050222001188385549510401465815187079080512838891
B=64019533680030876408443198762210829058751700634554282185987325820393598524794
S=K*pow(3,p-2,p)%p
O=None
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    x1,y1=P; x2,y2=Q
    if x1==x2:
        if (y1+y2)%p==0: return None
        l=3*x1*x1%p*pow(2*y1%p,p-2,p)%p
    else:
        l=(y2-y1)*pow((x2-x1)%p,p-2,p)%p
    x3=(l*l-x1-x2)%p; y3=(l*(x1-x3)-y1)%p
    return (x3,y3)
def neg(P): return None if P is None else (P[0],(-P[1])%p)
def mul(k,P):
    if k<0: return mul(-k,neg(P))
    R=None
    while k:
        if k&1: R=add(R,P)
        P=add(P,P); k>>=1
    return R
def oncurve(P):
    if P is None: return True
    x,y=P; return (y*y-pow(x,3,p)-B)%p==0
def load():
    BP=json.load(open('bitpoints.json'))
    P={int(b):((int(v[1])+S)%p,int(v[2])%p) for b,v in BP.items()}
    side={int(b):v[0] for b,v in BP.items()}
    return P,side
def target():
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    return ((C2+S)%p, C1%p)
def cornacchia_27(n):
    """solve 4n = L^2 + 27 M^2"""
    # find sqrt(-27) mod n ... use sqrt(-3): 4n = L^2+27M^2 <=> n = ((L+3M)/2)^2 + ... use standard
    # solve x^2 = -27 mod 4n via sqrt(-3) mod p
    import sympy
    return None
if __name__=='__main__':
    P,side=load(); T=target()
    print('all on curve:',all(oncurve(v) for v in P.values()),oncurve(T))
    print('B =',B)
    # group order candidates for j=0 curve: 4p = L^2+27M^2
    # find L,M
    import sympy
    # sqrt(-3) mod p
    r=sympy.ntheory.residue_ntheory.sqrt_mod(-3 % p, p)
    print('sqrt(-3) mod p ok:', (r*r+3)%p==0)
    # Cornacchia for x^2+3y^2 = p
    a0,b0=p,r%p
    lim=math.isqrt(p)
    while b0*b0>p:
        a0,b0=b0,a0%b0
    c=b0; d=math.isqrt((p-c*c)//3) if (p-c*c)%3==0 else None
    print('c=',c,'d=',d, 'check', None if d is None else c*c+3*d*d==p)
