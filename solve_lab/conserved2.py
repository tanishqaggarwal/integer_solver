#!/usr/bin/env python3
"""Extract the single conserved functional c (16-vector), its value, the 16 support inputs,
and whether the 2 MUX bits or any discrete flip can move it."""
import heal_harness as H
from jac_lib import D
import flint, json
p=H.p
def setstate(a=1,b=0):
    vA=H.loadd('best_agentA_39022.json')
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    H.val[2081]=a; H.val[4287]=b
    H.val[17325]=0; H.val[9413]=0; H.forward()
    H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
setstate()
leaf=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
freelist=sorted(H.freeinp); colidx={j:k for k,j in enumerate(freelist)}
def grads():
    vd=[None]*H.NVARS
    for j in H.freeinp: vd[j]=D(H.val[j],{colidx[j]:1})
    ns={'v':vd,'__builtins__':{}}
    for k,t in enumerate(H.order):
        r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
    R=[];G=[]
    for i in leaf:
        rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
        if isinstance(rr,D): R.append(rr.v%p); G.append(rr.g)
        else: R.append(rr%p); G.append({})
    return R,G
R,G=grads()
supp=sorted(set().union(*[set(g) for g in G]))
print('16 support free inputs (var ids):',[freelist[c] for c in supp])
# classify support: bits vs words vs leaves/slacks
bits=set(json.load(open('all_bits.json')))
pins=json.load(open('pinrec.json')); pintgt=set(pr[2] for pr in pins)
for c in supp:
    v=freelist[c]
    kind='BIT' if v in bits else ('PIN-TARGET' if v in pintgt else 'free')
    print(f'   x_{v}: {kind}')
# left null space of J (16 x |supp|): find c with c^T J=0
sidx={c:k for k,c in enumerate(supp)}
ctx=flint.fmpz_mod_ctx(p)
# Build J^T ( |supp| x 16 ), null space gives c
JT=flint.fmpz_mod_mat(len(supp),16,ctx)
for r,g in enumerate(G):
    for cc,co in g.items(): JT[sidx[cc],r]=co%p
# null space of J (right null of J is delta; we want left null = null of J^T? )
# left null of J: vectors c (16) with c^T J =0  <=> J^T c =0 => c in nullspace of J^T
null=JT.nullspace()[0]  # columns are basis of null(J^T) = left-null(J)
nd=null.ncols()
print('left-null(J) dim:',nd)
# find the combo among null basis with c^T R != 0
for jc in range(nd):
    c=[int(null[r,jc]) for r in range(16)]
    val=sum(c[r]*R[r] for r in range(16))%p
    tag='<== OBSTRUCTION' if val!=0 else '(annihilates R)'
    print(f'  null vec {jc}: c^T R = {val}  {tag}')
    if val!=0:
        nz=[(leaf[r],c[r]) for r in range(16) if c[r]]
        print('    combo (eq,coef):',nz)
