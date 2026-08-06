import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
atoms=[]; ateqs={}; reprs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        d=json.loads(line); atoms.append([(tuple(m),c) for m,c in d['poly']]); ateqs[i]=set(d.get('eqs',[])); reprs[i]=d.get('repr','')
def ev(poly,v):
    s=0
    for m,c in poly:
        t=c
        for var in m: t*=v[var]
        s+=t
    return s
# residue-shifted config
r8=109044024666698959972204451600908701898659086097062528124234304603594878834481
r9=33371159155735472537534252650716501592825364489306217536352743247010353604716
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[4287]=1
H.val[31861]=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
H.val[14865]=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
H.val[8731]=r8; H.val[9118]=r9; H.val[9413]=0; H.val[17325]=0
H.forward()
H.val[4432]=H.val[19964]; H.val[7068]=H.val[2099]
H.forward()
H.val[950]=H.val[9106]//(13523997*p); H.val[6947]=(6122989*H.val[2239])//p; H.val[33168]=-(H.val[31731]//p)
H.forward()
# examine 5 inconsistent rows
for eq in [24721,28737,29638,29959,37431]:
    print(f'=== eq {eq} ===')
    for ai,poly in enumerate(atoms):
        if eq in ateqs[ai]:
            v=ev(poly,H.val)
            if v!=0:
                # free inputs in this atom
                vs=set()
                for m,c in poly: vs.update(m)
                fr=sorted(vs & H.freeinp)
                print(f'  atom {ai} val%p={v%p!=0} :: {reprs[ai][:65]} | free:{fr}')
