import os, sys, random, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
SEVEN=[22229,22230,35758,35759,35760,35761,35762]
E=[2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125]
ES=set(E)
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(v,rounds=2):
    for _ in range(rounds):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,v)
            if nv is not None: v[u]=nv
    return v
v0=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
w=list(v0); fwd2(w,8)
avb=L.all_atom_values(w)
Ab=[avb[a] for a in SEVEN]

def setA(base, tgt):
    """Realise target A vector exactly, using the frame-2 handles.
       A0 = x_7068 - x_2099 - 7376877*x_642      (x_7068, x_642 free)
       A1 = x_28730 - p*x_9413                   (x_28730, x_9413 free)
       A2 = x_29854 - p*x_1329
       A3 = 5113045*x_9118 - x_29854
       A4 = x_31864 - p*x_10903
       A5 = x_8731 + x_31864
       A6 = x_642 - p*x_17325
    """
    v=list(base)
    A=tgt
    # --- A2, A3 : choose x_1329 and x_9118 so that A2+A3 = 5113045*x_9118 - p*x_1329
    S=A[2]+A[3]
    g=pow(5113045,-1,P)      # solve 5113045*x9118 == S (mod p) -> x1329 from remainder
    x9118=(S*g)%P
    rem=5113045*x9118 - S
    assert rem % P == 0
    x1329 = rem//P
    v[9118]=x9118; v[1329]=x1329
    v[29854]=A[2]+P*x1329
    # --- A4, A5
    D=A[5]-A[4]              # = x_8731 + p*x_10903
    x10903 = D//P
    x8731 = D - P*x10903
    v[10903]=x10903; v[8731]=x8731
    v[31864]=A[4]+P*x10903
    # --- A1 : x_28730 fixed mod p by base; choose x_9413
    #     A1 = x_28730 - p*x_9413 ; keep x_28730 as-is -> x_9413 = (x_28730-A1)/p
    num = base[28730]-A[1]
    if num % P: return None,'A1 not congruent to x_28730 mod p'
    v[9413]=num//P
    v[28730]=base[28730]
    # --- A6, A0 : A6 = x_642 - p*x_17325 ; A0 = x_7068 - x_2099 - 7376877*x_642
    #     keep x_7068 as-is; then x_642 = (x_7068 - x_2099 - A0)/7376877
    num2 = base[7068]-base[2099]-A[0]
    if num2 % 7376877: return None,'A0 not congruent mod 7376877'
    x642 = num2//7376877
    v[642]=x642; v[7068]=base[7068]
    num3 = x642 - A[6]
    if num3 % P: return None,'A6 mismatch mod p (congruence 1 violated)'
    v[17325]=num3//P
    fwd2(v,3)
    return v,None

# sanity: reproduce the delivered A exactly
v1,err = setA(w, Ab)
print('reproduce delivered A:', err)
if v1:
    av=L.all_atom_values(v1); f=set(L.failing_eqs(av))
    print('  score', L.NEQ-len(f), 'A matches:', [av[a] for a in SEVEN]==Ab, 'outside failures', sorted(f-ES))

# --- random A2..A5 with A0,A1,A6 kept: verify freedom ---
random.seed(7)
print('\n=== random A2..A5, A0/A1/A6 fixed ===')
for trial in range(6):
    tgt=list(Ab)
    for k in (2,3,4,5):
        tgt[k]=random.getrandbits(300)-(1<<299)
    vv,err=setA(w,tgt)
    if err: print(' err',err); continue
    av=L.all_atom_values(vv); f=set(L.failing_eqs(av))
    got=[av[a] for a in SEVEN]
    print(f' trial{trial}: A hit exactly={got==tgt}  score={L.NEQ-len(f)}  outside12_failures={sorted(f-ES)}')
