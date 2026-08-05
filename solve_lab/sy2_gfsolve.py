import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
fc=H.loadd('fc_partial.json')
# composite knob moves (maintain G1,G2 by slaving x_4432,x_7068)
COMPOSITE={8731:{8731:1,4432:1}, 9118:{9118:1,7068:1}}
FIXED={4287,2081,24601,31861,14865,4432,7068}  # don't move directly (pins & slaves)
def setup():
    for v in H.freeinp: H.val[v]=fc.get(v,0)
    H.forward()
setup()
F=set(H.fails())
# Build closure: iterate BFS eqs<->free inputs to fixpoint (bounded)
clo=set(F)
for _ in range(4):
    fr=set()
    for i in clo: fr|=(H.eqvars[i]&H.freeinp)
    fr-=FIXED
    new=set()
    for kn in fr:
        for i,vs in enumerate(H.eqvars):
            if kn in vs: new.add(i)
    if new<=clo: break
    clo|=new
    if len(clo)>1500: break
Feqs=sorted(clo)
knobs=set()
for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
knobs-=FIXED
knobs=sorted(knobs)
print('closure eqs:',len(Feqs),' knobs:',len(knobs))
def resids():
    ns={'v':H.val,'__builtins__':{}}
    return [eval(H.eqcode[i],ns)%p for i in Feqs]
base=resids()
nzidx=[k for k in range(len(Feqs)) if base[k]!=0]
print('nonzero rows:',len(nzidx))
def apply(mv,s=1):
    for k,v in mv.items(): H.val[k]=H.val[k]+s*v
def perturb_col(kn):
    mv=COMPOSITE.get(kn,{kn:1})
    apply(mv,1); H.forward(); Rn=resids(); apply(mv,-1); H.forward()
    return {k:(Rn[k]-base[k])%p for k in nzidx if (Rn[k]-base[k])%p!=0}
# build sparse columns
cols={}
for kn in knobs:
    c=perturb_col(kn)
    if c: cols[kn]=c
knobs=[k for k in cols]
print('effective knobs:',len(knobs))
# Clean GF(p) Gaussian: rows = nzidx, augmented with -base
# Represent rows as dict col->val plus rhs
rowdata={k:{} for k in nzidx}
for kn in knobs:
    for k,v in cols[kn].items(): rowdata[k][kn]=v
rhs={k:(-base[k])%p for k in nzidx}
rows=nzidx[:]
pivcol={}; pivrow_order=[]
free_rows=rows[:]
colset=knobs[:]
usedrows=set()
for kn in colset:
    prow=None
    for r in rows:
        if r in usedrows: continue
        if rowdata[r].get(kn,0)%p!=0: prow=r; break
    if prow is None: continue
    usedrows.add(prow); pivcol[kn]=prow; pivrow_order.append((kn,prow))
    ipv=inv(rowdata[prow][kn])
    # normalize pivot row
    for c in list(rowdata[prow]): rowdata[prow][c]=(rowdata[prow][c]*ipv)%p
    rhs[prow]=(rhs[prow]*ipv)%p
    # eliminate from all other rows
    for r in rows:
        if r==prow: continue
        f=rowdata[r].get(kn,0)%p
        if f==0: continue
        for c,val in rowdata[prow].items():
            rowdata[r][c]=(rowdata[r].get(c,0)-f*val)%p
        rhs[r]=(rhs[r]-f*rhs[prow])%p
incon=[r for r in rows if r not in usedrows and rhs[r]%p!=0 and all(v%p==0 for v in rowdata[r].values())]
print('pivots:',len(pivcol),' rank:',len(pivcol),' INCONSISTENT rows:',len(incon))
if incon:
    print('inconsistent eq indices:',[Feqs[r] for r in incon][:15])
else:
    print('*** mod-p system CONSISTENT ***')
