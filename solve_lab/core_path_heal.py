import pickle, json, time
import heal_harness as H
from collections import defaultdict
p=H.p
# load cone free inputs (must stay fixed to keep loads=0)
gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
def cone(roots):
    seen=set(); st=list(roots)
    while st:
        v=st.pop()
        if v in seen: continue
        seen.add(v)
        if v in gdef:
            for u in gdef[v][1]: st.append(u)
    return seen
loadcone=cone([11150,25739,37758,24908,35389,6671])  # x_24908 & loads & S,T
loadcone_free=loadcone & H.freeinp
print(f"load-cone free inputs to PIN: {len(loadcone_free)}")
# minimal core fix
d=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward(); V=H.val
V[14853]=V[12186]
V[16742]=(V[16742]//p)*p + V[24908]%p
H.forward()
F0=set(H.fails())
print(f"after minimal core fix: {len(F0)} fail, loads%p={[V[n]%p for n in [11150,25739,37758]]}")
# Now solve linear system for the ripple, PINNING loadcone_free to current residues.
# reuse lin_rels, add pins for loadcone_free (residue=current), add core already applied (x_14853,x_16742 pinned too)
lin_rels=pickle.load(open('lin_rels.pkl','rb'))
N=H.NVARS; CONST=N
par=list(range(N+1)); mult=[1]*(N+1); off=[0]*(N+1)
import sys; sys.setrecursionlimit(2000000)
def find2(x):
    ch=[]; r=x
    while par[r]!=r: ch.append(r); r=par[r]
    m,o=1,0
    for v in reversed(ch): m=(mult[v]*m)%p; o=(mult[v]*o+off[v])%p
    return r,m,o
contra=[]
def pin(v,k):
    rv,mv,ov=find2(v)
    if rv==CONST:
        if (mv+ov-k)%p: contra.append(v)
        return
    par[rv]=CONST; mult[rv]=0; off[rv]=((k-ov)*pow(mv,-1,p))%p
def union(a,b,m,o):
    ra,ma,oa=find2(a); rb,mb,ob=find2(b)
    if ra==rb:
        if (ma-m*mb)%p or (oa-(m*ob+o))%p: contra.append((a,b))
        return
    if ra==CONST: pin(b,(((ma+oa)-o)*pow(m,-1,p))%p); return
    if rb==CONST: pin(a,(m*(mb+ob)+o)%p); return
    inv=pow((m*mb)%p,-1,p); par[rb]=ra; mult[rb]=(ma*inv)%p; off[rb]=((oa-m*ob-o)*inv)%p
multi=[]
for ai,rel,const in lin_rels:
    vs=list(rel)
    if len(vs)==1: pin(vs[0],(-const*pow(rel[vs[0]],-1,p))%p)
    elif len(vs)==2:
        va,vb=vs; union(va,vb,(-rel[vb]*pow(rel[va],-1,p))%p,(-const*pow(rel[va],-1,p))%p)
    else: multi.append((rel,const))
# pin load-cone free inputs + core vars to current residues
for v in loadcone_free | {14853,12186,16742}:
    pin(v, V[v]%p)
print(f"unions+pins done, contra={len(contra)}")
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
    elif cc%p: contra.append('m')
print(f"reduced rels: {len(red)}, contra={len(contra)}")
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
print(f"gauss: pivots={len(pivots)}, contradictions={c2}")
if c2==0:
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
        if fv in loadcone_free or fv in {14853,12186,16742}: continue
        r,m,o=find2(fv)
        if r==CONST: nr=(m+o)%p
        elif r in rootval: nr=(m*rootval[r]+o)%p
        else: continue
        old=V[fv]; V[fv]=(old//p)*p+nr; ch+=(V[fv]%p!=old%p)
    print(f"changed {ch} non-loadcone free residues")
    H.forward()
    F=set(H.fails()); ns={'v':V,'__builtins__':{}}
    modp=sum(1 for i in F if eval(H.eqcode[i],ns)%p)
    print(f"AFTER ripple heal: {len(H.lines)-len(F)}/{len(H.lines)} sat, {len(F)} fail ({modp} mod-p)")
    print(f"loads%p={[V[n]%p for n in [11150,25739,37758]]}")
    json.dump({f'x_{i}':V[i] for i in range(N)},open('core_path_solved.json','w'))
