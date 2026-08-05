import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setfree(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
RIPPLE16=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
def rip(): return [i for i in RIPPLE16 if eval(H.eqcode[i],ns)!=0]
# step1: close gaps
setfree(7068,H.val[2099]); setfree(4432,H.val[19964])
print("after gap close, ripple fails:",len(rip()),rip())
# step2: close 7450 via x_2964, 7452 via x_24548
setfree(2964, H.val[26756]+H.val[579])
# 7452: 9367949*(x_24548-x_25442)-x_7927=0 -> x_24548 = x_25442 + x_7927*inv(9367949) but need exact Z. 
# x_7927 = p*x_11052; at base x_11052=? set x_24548 so 7452=0 in Z: 9367949*(x_24548-x_25442)=x_7927
inv=pow(9367949,-1,p)
# work mod p first: x_24548 = x_25442 + x_7927*inv (mod p) then adjust. Just set mod-p:
x7927=H.val[7927]; x25442=H.val[25442]
setfree(24548, (x25442 + x7927*inv)%p)
print("after closing 7450&7452:")
import json as J
for ai,expr in [('7450','v[2964]-v[26756]-v[579]'),('7452','9367949*(v[24548]-v[25442])-v[7927]'),
                ('44342_val',None),('45677_val',None)]:
    pass
atoms=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f: atoms.append(J.loads(line))
def ev(ai):
    s=0
    for mono,c in atoms[ai]['poly']:
        t=c
        for v in mono: t*=H.val[v]
        s+=t
    return s
for ai in [7450,7452,44342,45677]:
    print(f"  atom {ai} mod p = {ev(ai)%p}")
print("ripple fails now:",len(rip()),rip())
# Now check: is 44342/45677 a combo of 7450,7452? compute over 3 random continuous perturbations
import random
random.seed(1)
knobs=[2964,24548,26756,25442,7927,579]  # can't set gates; perturb free ancestors instead
print("\n-- test if 44342,45677 in span{7450,7452} via random free-input perturbations --")
freeperturb=[2964,24548,19569,11052,7068,4432]
def vec():
    return (ev(7450)%p, ev(7452)%p, ev(44342)%p, ev(45677)%p)
rows=[]
sv={v:H.val[v] for v in freeperturb}
for _ in range(6):
    for v in freeperturb: setfree(v, sv[v]+random.randint(-5,5))
    rows.append(vec())
    for v in freeperturb: setfree(v, sv[v])
import flint
ctx=flint.fmpz_mod_ctx(p)
M=flint.fmpz_mod_mat(len(rows),4,ctx)
for i,r in enumerate(rows):
    for j in range(4): M[i,j]=r[j]
print("rank of [7450,7452,44342,45677] samples:",M.rank(),"(if 2 => verifiers are combos)")
