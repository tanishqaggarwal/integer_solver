"""Cross-validate the s9/s10 atom model against the official checker.

Builds a state with random A2..A5 (should score 39,022 per lib) and writes it out
so `python3 checker.py` can be run on it independently.
"""
import os, sys, random, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E12=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds=3):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8)
avb=L.all_atom_values(w); Ab=[avb[a] for a in SEVEN]
def setA(base,A):
    v=list(base)
    S=A[2]+A[3]
    x9118=(S*pow(5113045,-1,P))%P
    rem=5113045*x9118-S
    assert rem%P==0
    x1329=rem//P
    v[9118]=x9118; v[1329]=x1329; v[29854]=A[2]+P*x1329
    D=A[5]-A[4]; x10903=D//P; x8731=D-P*x10903
    v[10903]=x10903; v[8731]=x8731; v[31864]=A[4]+P*x10903
    num=base[28730]-A[1]
    if num%P: return None
    v[9413]=num//P; v[28730]=base[28730]
    num2=base[7068]-base[2099]-A[0]
    if num2%7376877: return None
    x642=num2//7376877; v[642]=x642; v[7068]=base[7068]
    num3=x642-A[6]
    if num3%P: return None
    v[17325]=num3//P
    fwd2(v,3); return v
random.seed(11)
tgt=list(Ab)
for k in (2,3,4,5): tgt[k]=random.getrandbits(400)-(1<<399)
vv=setA(w,tgt)
av=L.all_atom_values(vv); f=L.failing_eqs(av)
print('lib says: score', L.NEQ-len(f), 'failing', sorted(f))
print('A realised exactly:', [av[a] for a in SEVEN]==tgt)
T.save(vv, os.path.join(HERE,'au_crossval.json'))
print('wrote au_crossval.json')
