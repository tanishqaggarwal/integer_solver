import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
def inv(x): return pow(x%p,-1,p)
fc=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0; H.val[950]=0; H.val[6947]=0
H.forward()
if H.val[37720]%9994531==0: H.val[8976]=H.val[37720]//9994531
H.forward()
# 2-hop closure of the 5 irreducible eqs, gather ALL fine-grained knobs (incl x_14865,x_31861 allowed)
seed=[29926,30383,32138,37297,38740]
clo=set(seed)
for _ in range(2):
    fr=set()
    for i in clo: fr|=(H.eqvars[i]&H.freeinp)
    fr-={4287,2081,24601,8731,9118,4432,7068,6418,12553}  # allow x_14865,x_31861 as knobs
    new=set()
    for kn in fr:
        for i,vs in enumerate(H.eqvars):
            if kn in vs: new.add(i)
    clo|=new
Feqs=sorted(clo)
knobs=set()
for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
knobs-={4287,2081,24601,8731,9118,4432,7068,6418,12553}
knobs=sorted(knobs)
def rc(idxs):
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns)%p for i in idxs}
base=rc(seed)
# only need to zero the 5 seed eqs; find fine knobs affecting them
cols={}
for kn in knobs:
    aff=[i for i in seed if kn in H.eqvars[i]]
    if not aff: continue
    H.val[kn]+=1; H.forward(); Rn=rc(seed); H.val[kn]-=1; H.forward()
    col={i:(Rn[i]-base[i])%p for i in seed if (Rn[i]-base[i])%p!=0}
    if col: cols[kn]=col
print('fine knobs touching the 5 irreducible eqs:',sorted(cols.keys()))
# Gaussian over the 5
rowdata={i:dict((kn,cols[kn][i]) for kn in cols if i in cols[kn]) for i in seed}
rhs={i:(-base[i])%p for i in seed}
used=set(); piv=0
for kn in cols:
    prow=None
    for i in seed:
        if i in used: continue
        if rowdata[i].get(kn,0)%p!=0: prow=i;break
    if prow is None: continue
    used.add(prow); piv+=1
    ipv=inv(rowdata[prow][kn])
    for c in list(rowdata[prow]): rowdata[prow][c]=(rowdata[prow][c]*ipv)%p
    rhs[prow]=(rhs[prow]*ipv)%p
    for i in seed:
        if i==prow: continue
        f=rowdata[i].get(kn,0)%p
        if f==0: continue
        for c,val in rowdata[prow].items(): rowdata[i][c]=(rowdata[i].get(c,0)-f*val)%p
        rhs[i]=(rhs[i]-f*rhs[prow])%p
incon=[i for i in seed if i not in used and rhs[i]%p!=0]
print('over 5 irreducible eqs: rank=%d inconsistent=%d %s'%(piv,len(incon),incon))
# save the loads-zeroed config as artifact
out={f'x_{i}':H.val[i] for i in range(H.NVARS) if H.val[i]!=0}
json.dump(out,open('sy2_loadszeroed_39018.json','w'))
print('saved sy2_loadszeroed_39018.json, fails=',len(H.fails()))
