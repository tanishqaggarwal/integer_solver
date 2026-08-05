import pickle, json, time
import heal_harness as H
from collections import defaultdict
p=H.p
R=pickle.load(open('lin_reduced.pkl','rb'))
par,mult,off,red_fixed=R['par'],R['mult'],R['off'],R['red']
CONST=H.NVARS
def find2(x):
    chain=[]; r=x
    while par[r]!=r: chain.append(r); r=par[r]
    m,o=1,0
    for v in reversed(chain): m=(mult[v]*m)%p; o=(mult[v]*o+off[v])%p
    return r,m,o
# wire set
import heal_harness
# genuine value*value gates: parse gate defs, product of two non-wire vars
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
# determine wire members: value %p==0 persistently. Use: a var is wire if its find2 root maps it with the 220-set. Reload wire from lin_extract logic quickly:
# Simpler: wire = vars whose forward value %p==0 AND appear as pure p. Use gates.jsonl constant gate x_26064=p + identity closure done earlier -> recompute:
par2=list(range(H.NVARS)); 
def f3(x):
    while par2[x]!=x: par2[x]=par2[par2[x]]; x=par2[x]
    return x
atoms=[json.loads(l) for l in open('atoms/poly_atoms.jsonl')]
for a in atoms:
    poly=a['poly']
    if all(len(vs)<=1 for vs,c in poly):
        lin=[(vs[0],c) for vs,c in poly if len(vs)==1]; const=sum(c for vs,c in poly if len(vs)==0)
        if len(lin)==2 and const==0 and abs(lin[0][1])==1 and abs(lin[1][1])==1:
            ra,rb=f3(lin[0][0]),f3(lin[1][0]); 
            if ra!=rb: par2[rb]=ra
wire=set(v for v in range(H.NVARS) if f3(v)==f3(26064))
# value*value gates (both factors non-wire)
vv_gates=[]
import re
VAR=re.compile(r'x_(\d+)')
for t,(rhs,vids) in gdef.items():
    if '*' in rhs:
        vs=[int(m) for m in VAR.findall(rhs)]
        # get the two product factors: gate is single product coef*(x_a*x_b) or x_a*x_b
        # factors = the two vars in the product (crude but gates are simple products)
        nonwire=[v for v in vs if v not in wire]
        if len(vs)==2 and vs[0] not in wire and vs[1] not in wire:
            vv_gates.append((t,vs[0],vs[1]))
print(f"genuine value*value gates: {len(vv_gates)}")

def reduce_row(pairs, k):
    r=defaultdict(int); cc=k%p
    for v,coef in pairs:
        rv,m,o=find2(v)
        if rv==CONST: cc=(cc+coef*(m+o))%p
        else: r[rv]=(r[rv]+coef*m)%p; cc=(cc+coef*o)%p
    return {kk:vv%p for kk,vv in r.items() if vv%p}, cc%p

def gauss(rows,consts):
    col_rows=defaultdict(set)
    for i,row in enumerate(rows):
        for v in row: col_rows[v].add(i)
    alive=set(range(len(rows))); pivots={}; contra=0
    for i in sorted(range(len(rows)), key=lambda i:len(rows[i])):
        if i not in alive: continue
        row=rows[i]
        if not row:
            if consts[i]%p: contra+=1
            alive.discard(i); continue
        pc=min(row, key=lambda v: len(col_rows[v]&alive))
        inv=pow(row[pc],-1,p); pivots[pc]=i; alive.discard(i)
        rows[i]={v:(c*inv)%p for v,c in row.items()}; consts[i]=(consts[i]*inv)%p; row=rows[i]
        for j in list(col_rows[pc]):
            if j==i or j not in alive: continue
            fj=rows[j].get(pc,0)
            if not fj: continue
            for v,c in row.items():
                nv=(rows[j].get(v,0)-fj*c)%p
                if nv: rows[j][v]=nv; col_rows[v].add(j)
                elif v in rows[j]: del rows[j][v]; col_rows[v].discard(j)
            consts[j]=(consts[j]-fj*consts[i])%p
            if not rows[j] and consts[j]%p: contra+=1; alive.discard(j)
    return pivots,contra

# load 39013
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
print(f"start fails: {len(H.fails())}")
for it in range(8):
    V=H.val
    rows=[dict(r) for r,c in red_fixed]; consts=[c for r,c in red_fixed]
    # value*value linearized: x_out - b0*x_a - a0*x_b + a0*b0 = 0
    for t,va,vb in vv_gates:
        a0=V[va]%p; b0=V[vb]%p
        r,cc=reduce_row([(t,1),(va,(-b0)%p),(vb,(-a0)%p)], (a0*b0)%p)
        if r: rows.append(r); consts.append(cc)
        elif cc%p: pass
    # core conditions
    for a,b in [(14853,12186),(24908,16742)]:
        r,cc=reduce_row([(a,1),(b,-1)],0)
        if r: rows.append(r); consts.append(cc)
    pivots,contra=gauss(rows,consts)
    # back-substitute: free roots keep current residue
    allroots=set()
    for row in rows: allroots|=set(row)
    rootval={rt:(V[rt]%p if rt<H.NVARS else 1) for rt in allroots if rt not in pivots}
    for col,ri in sorted(pivots.items(), key=lambda kv:-kv[1]):
        row=rows[ri]; s=consts[ri]
        for v,c in row.items():
            if v!=col: s=(s-c*rootval.get(v,0))%p
        rootval[col]=s%p
    # apply to free inputs
    for fv in H.freeinp:
        r,m,o=find2(fv)
        if r==CONST: nr=(m+o)%p
        elif r in rootval: nr=(m*rootval[r]+o)%p
        else: continue
        old=V[fv]; q=old//p
        V[fv]=q*p+nr
    H.forward()
    F=H.fails()
    ns={'v':V,'__builtins__':{}}
    modp=sum(1 for i in F if eval(H.eqcode[i],ns)%p!=0)
    print(f"iter {it}: {len(H.lines)-len(F)}/{len(H.lines)} sat, {len(F)} fail ({modp} mod-p, {len(F)-modp} carry), contra={contra}",flush=True)
    if len(F)==0: json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('NEWTON_SOLVED.json','w')); print("*** SOLVED ***"); break
    if modp==0: print("all mod-p satisfied - only carries left!"); json.dump({f'x_{i}':V[i] for i in range(H.NVARS)},open('newton_modp_done.json','w')); break
