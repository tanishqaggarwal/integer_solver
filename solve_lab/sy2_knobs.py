import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
fc=H.loadd('fc_partial.json')
COMPOSITE={8731:{8731:1,4432:1}, 9118:{9118:1,7068:1}}
FIXED={4287,2081,24601,31861,14865,4432,7068}
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.forward()
F=sorted(H.fails())
print('fails:',len(F),F)
# closure (2 hops enough for gadgets)
clo=set(F)
for _ in range(3):
    fr=set()
    for i in clo: fr|=(H.eqvars[i]&H.freeinp)
    fr-=FIXED
    new=set()
    for kn in fr:
        for i,vs in enumerate(H.eqvars):
            if kn in vs: new.add(i)
    if new<=clo: break
    clo|=new
Feqs=sorted(clo)
print('closure eqs:',len(Feqs))
def resids():
    ns={'v':H.val,'__builtins__':{}}
    return [eval(H.eqcode[i],ns)%p for i in Feqs]
base=resids()
nz=[Feqs[k] for k in range(len(Feqs)) if base[k]!=0]
print('nonzero rows:',len(nz),nz)
knobs=set()
for i in Feqs: knobs|=(H.eqvars[i]&H.freeinp)
knobs-=FIXED
# find effective: perturb, check if changes any nonzero row (sparse eval)
eff=[]
nzset=set(nz)
for kn in sorted(knobs):
    mv=COMPOSITE.get(kn,{kn:1})
    # only eqs containing kn (or its composite parts) can change
    parts=set(mv)
    aff=[i for i in Feqs if H.eqvars[i]&parts]
    for k,v in mv.items(): H.val[k]+=v
    H.forward()
    ns={'v':H.val,'__builtins__':{}}
    changed=any(eval(H.eqcode[i],ns)%p != base[Feqs.index(i)] for i in aff)
    for k,v in mv.items(): H.val[k]-=v
    H.forward()
    if changed: eff.append(kn)
print('effective knobs:',len(eff),eff)
