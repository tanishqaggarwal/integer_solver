import os, sys, math, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v=L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
print('x_2099 == x_6418 :', v[2099]==v[6418])
print('x_19964 == x_12553:', v[19964]==v[12553])
print('gcd(7376877,p) =', math.gcd(7376877,P), ' gcd(5113045,p) =', math.gcd(5113045,P),
      ' gcd(15804267,p) =', math.gcd(15804267,P))
print('p prime?', pow(2,P-1,P)==1)
# a37887 is a perfect square of a linear form: check numerically on random perturbations
DETACH={7068:22229,28730:22230,29854:35758,31864:35761,642:35762}
definer={t:a for t,a in L.definer.items() if t not in DETACH}
ORDER=[t for t in ad.ORDER if t not in DETACH]
def fwd2(x,r=3):
    for _ in range(r):
        for u in ORDER:
            nv=T.solve_lin(definer[u],u,x)
            if nv is not None: x[u]=nv
    return x
w=list(v); fwd2(w,8)
import random
random.seed(5)
ok=0; bad=0
for t in range(8):
    x=list(w)
    for u in (4432,19964,28730,18253,9629,23754,35619,7945,23642,23822,37254,15324,37720,30108,34600):
        x[u]+= random.randint(-5,5)
    val=L.evalpoly(L.polys[37887],x)
    r=math.isqrt(abs(val))
    if val>=0 and r*r==val: ok+=1
    else: bad+=1
print(f'a37887 a perfect square on {ok}/{ok+bad} random perturbations (bad={bad})')
m,sq,co=L.eq_atoms[8680]
print('eq 8680: mult',m,'is_square',sq,'atoms',co)
# every equation that a22231 touches is in E12
E12=set([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
print('eqs(a22231) subset of E12:', set(L.atom2eq[22231])<=E12, sorted(L.atom2eq[22231]))
