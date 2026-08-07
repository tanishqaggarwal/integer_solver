"""Which upstream residues (S,T,U,B,A) make a small drop-set feasible?

local1 proved 7 is optimal with S,T,U,B,A frozen.  Those five numbers are the ONLY channel
through which the rest of the circuit speaks to this neighbourhood:

    S = x1956*x17065     T = 5113045*x7075*x9118     U = x7075*x8731
    A = x7068 - x2099    B = x28730

Feasibility of a kept-set K is  (S,T,U,B,A) in phi(ker M_K), where
phi(y) = (z0+z1, y3+y4, y6-y5, y2, y1+7376877*y7) mod (p,p,p,p,7376877p).
So for each K we can read off EXACTLY what the rest of the circuit must deliver.
"""
import sys, os, itertools, time
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw, atomval
P = L.P; Q = 7376877
HERE = os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
v = load_raw(os.path.join(LAB,'best','new_instance_partial_39026.json'))
AV=[atomval(a,v) for a in range(L.NA)]
KEY=[35756,35757,35758,35759,35760,35761,35762,22229,22230]
EQS=sorted(set().union(*[set(L.atom2eq.get(a,{})) for a in KEY]))
rows=[[L.eq_atoms[e][2].get(a,0) for a in KEY] for e in EQS]

def rank(M):
    M=[[Fraction(x) for x in r] for r in M]; r=0; n=len(M[0])
    for c in range(n):
        pv=next((i for i in range(r,len(M)) if M[i][c]),None)
        if pv is None: continue
        M[r],M[pv]=M[pv],M[r]
        pr=M[r]; inv=Fraction(1,1)/pr[c]
        M[r]=[x*inv for x in pr]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[M[i][k]-f*M[r][k] for k in range(n)]
        r+=1
    return r
print("rank of the full 15x9 system:", rank(rows), " (9 => the only local solution is all-atoms-zero)")
print()
print("=> A FULL local fix (0 failures here) requires, from the rest of the circuit:")
print("     S = x1956*x17065        == 0  (mod p)      currently", (v[1956]*v[17065])%P==0)
print("     T = 5113045*x7075*x9118 == 0  (mod p)  <=>  p | x9118      currently", v[9118]%P==0)
print("     U = x7075*x8731         == 0  (mod p)  <=>  p | x8731      currently", v[8731]%P==0)
print("     B = x28730              == 0  (mod p)                      currently", v[28730]%P==0)
print("     A = x7068 - x2099       == 0  (mod 7376877*p)              currently",
      (v[7068]-v[2099])%(Q*P)==0)
