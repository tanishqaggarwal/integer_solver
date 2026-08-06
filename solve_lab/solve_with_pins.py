import pickle, json, time
import heal_harness as H
from collections import defaultdict
p=H.p
# check if core vars are pinned
pinrec=json.load(open('pinrec.json'))
pintgts={tgt for i,sel,tgt,const,coef,handle in pinrec}
for v in [14853,12186,16742,24908]:
    print(f"x_{v} is pin target: {v in pintgts}")
# load 39013, get active selectors
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
# build extended linear relations: fixed reduced + active pins
lin_rels=pickle.load(open('lin_rels.pkl','rb'))
# add active pins as 1-term relations target ≡ CONST
pin_rels=[]
active=0
for i,sel,tgt,const,coef,handle in pinrec:
    s=V[sel]%p
    if s==1:
        pin_rels.append(({tgt:1}, (-const)%p)); active+=1   # x_tgt - CONST = 0 -> {tgt:1}, const=-CONST
    # inactive: handle partner ≡0 (skip - handles are separate)
print(f"active pins: {active}")
all_rels=[(None,rel,c) for (rel,c) in [(r,c) for r,c in [(dict(rr),cc) for rr,cc in [(x[0],x[1]) for x in pin_rels]]]]
# simpler: rebuild combined rel list
combined=[]
for ai,rel,const in lin_rels: combined.append((rel,const))
for rel,const in pin_rels: combined.append((rel,const))
# union-find + reduce (reuse logic)
N=H.NVARS; CONST=N
par=list(range(N+1)); mult=[1]*(N+1); off=[0]*(N+1)
import sys; sys.setrecursionlimit(2000000)
def find2(x):
    chain=[]; r=x
    while par[r]!=r: chain.append(r); r=par[r]
    m,o=1,0
    for v in reversed(chain): m=(mult[v]*m)%p; o=(mult[v]*o+off[v])%p
    return r,m,o
contra=[]
def pin(v,k):
    rv,mv,ov=find2(v)
    if rv==CONST:
        if (mv+ov-k)%p: contra.append(('pin',v))
        return
    par[rv]=CONST; mult[rv]=0; off[rv]=((k-ov)*pow(mv,-1,p))%p
def union(a,b,m,o):
    ra,ma,oa=find2(a); rb,mb,ob=find2(b)
    if ra==rb:
        if (ma-m*mb)%p or (oa-(m*ob+o))%p: contra.append(('u',a,b)); 
        return
    if ra==CONST: pin(b,(((ma+oa)-o)*pow(m,-1,p))%p); return
    if rb==CONST: pin(a,(m*(mb+ob)+o)%p); return
    inv=pow((m*mb)%p,-1,p); par[rb]=ra; mult[rb]=(ma*inv)%p; off[rb]=((oa-m*ob-o)*inv)%p
multi=[]
for rel,const in combined:
    vs=list(rel)
    if len(vs)==1: v=vs[0]; pin(v,(-const*pow(rel[v],-1,p))%p)
    elif len(vs)==2:
        va,vb=vs; union(va,vb,(-rel[vb]*pow(rel[va],-1,p))%p,(-const*pow(rel[va],-1,p))%p)
    else: multi.append((rel,const))
print(f"unions done, contra={len(contra)}, multi={len(multi)}")
def reduce_row(rel,const):
    r=defaultdict(int); cc=const%p
    for v,c in rel.items():
        rv,m,o=find2(v)
        if rv==CONST: cc=(cc+c*(m+o))%p
        else: r[rv]=(r[rv]+c*m)%p; cc=(cc+c*o)%p
    return {k:v%p for k,v in r.items() if v%p}, cc%p
red=[]
for rel,const in multi:
    r,cc=reduce_row(rel,const)
    if r: red.append((r,cc))
    elif cc%p: contra.append(('m',))
# core conditions
for a,b in [(14853,12186),(24908,16742)]:
    r,cc=reduce_row({a:1,b:-1},0)
    if r: red.append((r,cc))
    elif cc%p: contra.append(('core',a,b))
print(f"reduced+core rels: {len(red)}, contra={len(contra)}")
# gauss
rows=[dict(r) for r,c in red]; consts=[c for r,c in red]
col_rows=defaultdict(set)
for i,row in enumerate(rows):
    for v in row: col_rows[v].add(i)
alive=set(range(len(rows))); pivots={}; c2=0
for i in sorted(range(len(rows)),key=lambda i:len(rows[i])):
    if i not in alive: continue
    row=rows[i]
    if not row:
        if consts[i]%p: c2+=1
        alive.discard(i); continue
    pc=min(row,key=lambda v:len(col_rows[v]&alive)); inv=pow(row[pc],-1,p); pivots[pc]=i; alive.discard(i)
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
        if not rows[j] and consts[j]%p: c2+=1; alive.discard(j)
print(f"gauss pivots={len(pivots)}, contradictions={c2}")
if c2>0: print("INCONSISTENT with pins+core"); 
else:
    print("CONSISTENT with pins+core!")
    allroots=set()
    for row in rows: allroots|=set(row)
    rootval={rt:(V[rt]%p if rt<N else 1) for rt in allroots if rt not in pivots}
    for col,ri in sorted(pivots.items(),key=lambda kv:-kv[1]):
        row=rows[ri]; s=consts[ri]
        for v,c in row.items():
            if v!=col: s=(s-c*rootval.get(v,0))%p
        rootval[col]=s%p
    ch=0
    for fv in H.freeinp:
        r,m,o=find2(fv)
        if r==CONST: nr=(m+o)%p
        elif r in rootval: nr=(m*rootval[r]+o)%p
        else: continue
        old=V[fv]; V[fv]=(old//p)*p+nr; ch+= (V[fv]%p!=old%p)
    print(f"changed {ch} free residues")
    H.forward()
    F=H.fails(); ns={'v':V,'__builtins__':{}}
    modp=sum(1 for i in F if eval(H.eqcode[i],ns)%p)
    print(f"AFTER pins+core solve: {len(H.lines)-len(F)}/{len(H.lines)} sat, {len(F)} fail ({modp} mod-p)")
    print(f"x_29322%p={(V[14853]-V[12186])%p}, x_3558%p={(V[24908]-V[16742])%p}")
    json.dump({f'x_{i}':V[i] for i in range(N)},open('pins_core_solved.json','w'))
