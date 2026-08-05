import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
from collections import defaultdict
p=H.p
# load atoms
atoms=[]; ateqs=[]
with open('atoms/poly_atoms.jsonl') as f:
    for line in f:
        d=json.loads(line); atoms.append([(tuple(m),c) for m,c in d['poly']]); ateqs.append(set(d.get('eqs',[])))
# map free var -> atoms containing it (as a constraint, i.e. where var can be a knob)
FREE=H.freeinp
FIXED={4287,2081,24601,8731,9118,31861,14865,4432,7068}  # regime, load-residues, pins, slaved
def atom_val(poly):
    v=H.val; s=0
    for m,c in poly:
        t=c
        for var in m: t*=v[var]
        s+=t
    return s
def coeff_of(poly,var):
    # linear coefficient of var at current point (assumes var appears linearly)
    v=H.val; c=0
    for m,c0 in poly:
        cnt=m.count(var)
        if cnt==1:
            t=c0
            for var2 in m:
                if var2!=var: t*=v[var2]
            c+=t
        elif cnt>1:
            return None  # nonlinear in var
    return c
# precompute: for each free var, list of atom indices where it appears
var2at=defaultdict(list)
for ai,poly in enumerate(atoms):
    seen=set()
    for m,c in poly:
        for var in m:
            if var in FREE and var not in seen:
                seen.add(var); var2at[var].append(ai)
# healability: free var v can cancel atom ai if coeff c divides residual (exact) 
def setup():
    for v in FREE: H.val[v]=vA.get(v,0)
    H.val[4287]=1
    H.val[31861]=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
    H.val[14865]=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
    H.val[8731]=r8; H.val[9118]=r9; H.val[9413]=0; H.val[17325]=0
    H.forward(); slave()
def slave():
    H.val[4432]=H.val[19964]+H.val[28730]
    H.val[7068]=7376877*H.val[642]+H.val[2099]
    H.forward()
r8=109044024666698959972204451600908701898659086097062528124234304603594878834481
r9=33371159155735472537534252650716501592825364489306217536352743247010353604716
vA=H.loadd('best_agentA_39022.json')
setup()
# load handles each pass
def set_load_handles():
    if H.val[9106]%(13523997*p)==0 or True:
        H.val[950]=H.val[9106]//(13523997*p)
    H.val[6947]=(6122989*H.val[2239])//p
    H.val[33168]=-(H.val[31731]//p)
set_load_handles(); H.forward(); slave()
print('start fails:',len(H.fails()))
# iterative structural heal
for it in range(40):
    F=H.fails()
    if not F: print('SOLVED at iter',it); break
    # find nonzero atoms in failing eqs
    Fset=set(F)
    nz=[]
    for ai,poly in enumerate(atoms):
        if ateqs[ai] & Fset:
            if atom_val(poly)!=0: nz.append(ai)
    # heal each nonzero atom by a free knob (prefer knobs in fewest atoms, not fixed)
    healed_vars=set(); nhealed=0
    # sort atoms by number of free knobs (fewest options first)
    def knobs_for(ai):
        poly=atoms[ai]; V=atom_val(poly); res=[]
        for var in var2at_local(ai):
            if var in FIXED or var in healed_vars: continue
            c=coeff_of(poly,var)
            if c is None or c==0: continue
            if V % c==0: res.append((var,c))
        return res
    def var2at_local(ai):
        s=set()
        for m,c in atoms[ai]:
            for var in m:
                if var in FREE: s.add(var)
        return s
    cand=[]
    for ai in nz:
        ks=knobs_for(ai)
        if ks: cand.append((ai,ks))
    # greedily heal: pick knob appearing in fewest OTHER atoms
    cand.sort(key=lambda x: min(len(var2at[v]) for v,_ in x[1]))
    for ai,ks in cand:
        # choose knob with fewest atom memberships
        ks2=[(v,c) for v,c in ks if v not in healed_vars]
        if not ks2: continue
        ks2.sort(key=lambda vc: len(var2at[vc[0]]))
        var,c=ks2[0]
        V=atom_val(atoms[ai])
        if V%c!=0: continue
        H.val[var]=H.val[var]-V//c
        healed_vars.add(var); nhealed+=1
    set_load_handles(); H.forward(); slave()
    nf=len(H.fails())
    print(f'iter {it}: healed {nhealed} atoms -> fails={nf}')
    if nf==0: break
print('final fails:',len(H.fails()))
# save if improved
F=H.fails()
if len(F)<11:
    out={f'x_{i}':H.val[i] for i in range(H.NVARS) if H.val[i]!=0}
    json.dump(out,open('sy_healer_out.json','w'))
    print('SAVED sy_healer_out.json with',len(F),'fails')
