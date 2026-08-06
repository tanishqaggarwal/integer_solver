import pickle, sys, time
import heal_harness as H
from collections import defaultdict
p=H.p
R=pickle.load(open('lin_reduced.pkl','rb'))
par,mult,off=R['par'],R['mult'],R['off']
G=pickle.load(open('gauss_result.pkl','rb'))
rows,consts,pivots=G['rows'],G['consts'],G['pivots']
CONST=H.NVARS
def find2(x):
    chain=[]; r=x
    while par[r]!=r: chain.append(r); r=par[r]
    m,o=1,0
    for v in reversed(chain): m=(mult[v]*m)%p; o=(mult[v]*o+off[v])%p
    return r,m,o
# current residues at 39013
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
cur_res={}  # root -> current residue
def root_res(rt):
    # find a var mapping to this root to read current residue: use rt itself
    return H.val[rt]%p if rt<H.NVARS else 1
# back-substitution: free roots -> current residue; pivots in reverse order
pivot_cols=set(pivots)
all_roots=set()
for row in rows: all_roots|=set(row)
free_roots=[rt for rt in all_roots if rt not in pivot_cols]
rootval={}
for rt in free_roots: rootval[rt]=H.val[rt]%p if rt<H.NVARS else 1
# order pivots by elimination order (row index) then reverse
piv_by_row=sorted(pivots.items(), key=lambda kv: kv[1])  # (col,rowidx)
for col,ri in reversed(piv_by_row):
    row=rows[ri]; s=consts[ri]
    for v,coef in row.items():
        if v==col: continue
        s=(s-coef*rootval.get(v,0))%p
    rootval[col]=(s*pow(row[col],-1,p))%p if row[col]!=1 else s%p
print(f"solved {len(rootval)} root residues")
# apply to free inputs: set residue, keep quotient
changed=0
for fv in H.freeinp:
    r,m,o=find2(fv)
    if r==CONST: newres=(m+o)%p
    elif r in rootval: newres=(m*rootval[r]+o)%p
    else: continue  # not in system, leave
    old=H.val[fv]; q=old//p if old>=0 else -((-old+p-1)//p)
    H.val[fv]=q*p+newres
    if H.val[fv]%p!=old%p: changed+=1
print(f"changed {changed} free-input residues")
H.forward()
F=H.fails()
print(f"AFTER LINEAR SOLVE: {len(H.lines)-len(F)}/{len(H.lines)} satisfied ({len(F)} fail)")
# check mod-p vs carry among fails
ns={'v':H.val,'__builtins__':{}}
modp=sum(1 for i in F if eval(H.eqcode[i],ns)%p!=0); carry=len(F)-modp
print(f"  of {len(F)} fails: {modp} fail mod p, {carry} fail only in Z (carry)")
# core check
print(f"x_29322%p={(H.val[14853]-H.val[12186])%p}, x_3558%p={(H.val[24908]-H.val[16742])%p}")
import json
json.dump({f'x_{i}':H.val[i] for i in range(H.NVARS)}, open('linear_solved.json','w'))
print("saved linear_solved.json")
