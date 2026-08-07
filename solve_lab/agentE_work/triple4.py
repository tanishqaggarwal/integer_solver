"""Exact two-unknown congruence solve of the obstruction triple, using the two knobs on which
   (U,V) is exactly affine:  x_30468 and x_33169."""
import sys, json, math
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(3000000)
import engine as E, fast
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
U_IDX,V_IDX=29210,8736
K1,K2=30468,33169
s={int(k):int(v) for k,v in json.load(open('triple_state_seed.json')).items()}
v0=E.forward(s); U0,V0=v0[U_IDX],v0[V_IDX]
def col(f):
    o=v0[f]
    v1,_=fast.apply_delta(v0,{f:o+1}); v2,_=fast.apply_delta(v0,{f:o+2}); v3,_=fast.apply_delta(v0,{f:o+7})
    a1=v1[U_IDX]-U0; b1=v1[V_IDX]-V0
    lin = (v2[U_IDX]-U0==2*a1 and v2[V_IDX]-V0==2*b1 and v3[U_IDX]-U0==7*a1 and v3[V_IDX]-V0==7*b1)
    return a1,b1,lin
A,Bc,l1=col(K1); Cc,D,l2=col(K2)
print("x_%d linear over {1,2,7}: %s ; x_%d: %s"%(K1,l1,K2,l2))
print("dU/dK1 bits=%d  dV/dK1 bits=%d  dU/dK2=%s  dV/dK2 bits=%d"%(A.bit_length(),Bc.bit_length(),Cc,D.bit_length()))
assert Cc==0
M=5002401*A+15322661*Bc
N0=5002401*U0+15322661*V0
MOD=15322661*D
print("gcd(A,p) =",math.gcd(A,P))
# condition 1:  U0 + A*d1 = 0 (mod p)
g=math.gcd(A,P)
if U0%g:
    print("CONDITION 1 UNSOLVABLE"); sys.exit()
d1_0=(-U0*pow(A//g,-1,P//g))%(P//g)
step=P//g
print("d1 = %d + k*%d   (bits %d)"%(d1_0%step,step,d1_0.bit_length()))
# condition 2:  N0 + M*d1 + 15322661*D*d2 = 0  -> need MOD | (N0 + M*d1)
# with d1 = d1_0 + k*step :  N0 + M*d1_0 + M*step*k = 0 (mod MOD)
a=(M*step)%MOD; b=(-(N0+M*d1_0))%MOD
g2=math.gcd(a,MOD)
print("condition 2: %d*k = %d (mod %d);  gcd=%d ; divides rhs? %s"%(a and 1,b and 1,MOD and 1,g2,b%g2==0))
if b%g2:
    print("NO SOLUTION on this 2-knob subspace"); sys.exit()
k=(b//g2)*pow(a//g2,-1,MOD//g2)%(MOD//g2)
d1=d1_0+k*step
num=N0+M*d1
assert num%MOD==0
d2=-num//MOD
print("SOLUTION FOUND: d1 bits=%d  d2 bits=%d"%(d1.bit_length(),d2.bit_length()))
ns=dict(s); ns[K1]=v0[K1]+d1; ns[K2]=v0[K2]+d2
v=E.forward(ns); U,V=v[U_IDX],v[V_IDX]
print("check: p|U ->",U%P==0,"  p|V ->",V%P==0,"  exact row 5002401U+15322661V==0 ->",5002401*U+15322661*V==0)
av=E.badatoms(v); ff=E.eqfails(av)
print("EXACT: fails=%d score=%d  bad=%s"%(len(ff),39033-len(ff),sorted(av)[:20]))
json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open('triple4_%d.json'%(39033-len(ff)),'w'))
json.dump({str(a2):str(int(b2)) for a2,b2 in ns.items()}, open('triple4_seed.json','w'))
