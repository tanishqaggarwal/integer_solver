"""Complete the x4287 route: also get 13523997*p | x9106.

Shifting x9118 by p*j and x8731 by p*k leaves x27177 and x4306 zero mod p (their responses are
affine, so the shift moves them by multiples of p), but it does move x9106 -- affinely.  So fit
x9106/p as an affine function of (j,k) and solve the single congruence mod 13523997.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; M=13523997
HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
base=load_raw(os.path.join(HERE,'data','sw6_out.json'))
n9, n8 = base[9118], base[8731]
def shift(j,k):
    v=list(base); L.ripple(v,{9118:n9+P*j, 8731:n8+P*k}); return v
def probe(j,k):
    v=shift(j,k)
    assert v[27177]%P==0 and v[4306]%P==0, "mod-p work disturbed"
    assert v[9106]%P==0
    return (v[9106]//P) % M
f00=probe(0,0); f10=probe(1,0); f01=probe(0,1)
a=(f10-f00)%M; b=(f01-f00)%M
print(f"x9106/p mod {M}:  f(0,0)={f00}  d/dj={a}  d/dk={b}")
import math
g=math.gcd(math.gcd(a,b),M)
print(f"gcd(a,b,M)={g}   solvable: {f00%g==0}")
sol=None
if f00%g==0:
    # solve a*j + b*k == -f00 (mod M)
    for k in range(0, M//math.gcd(b,M) if b else 1):
        r=(-f00-b*k)%M
        d=math.gcd(a,M)
        if r%d==0:
            j=(r//d)*pow(a//d,-1,M//d)%(M//d)
            sol=(j,k); break
        if k>200000: break
print("solution (j,k):",sol)
if sol:
    j,k=sol
    v=shift(j,k)
    print("  13523997*p | x9106 :", v[9106]%(M*P)==0)
    L.ripple(v,{950: v[9106]//(M*P), 6947: 6122989*v[2239]//P, 33168: -v[31731]//P})
    AV=[L.evalpoly(L.polys[a2],v) for a2 in range(L.NA)]
    F=L.failing_eqs(AV); B=[a2 for a2 in range(L.NA) if AV[a2]!=0]
    print(f"  RESULT broken atoms={B} failing={len(F)} score={L.NEQ-len(F)}")
    json.dump({('x_%d'%i):v[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','sw7_out.json'),'w'))
