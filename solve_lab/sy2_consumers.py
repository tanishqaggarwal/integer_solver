import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
atoms=[]; ateqs=[]; reprs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        d=json.loads(line); atoms.append([(tuple(m),c) for m,c in d['poly']]); ateqs.append(set(d.get('eqs',[]))); reprs[i]=d.get('repr','')
def ev(poly):
    v=H.val; s=0
    for m,c in poly:
        t=c
        for var in m: t*=v[var]
        s+=t
    return s
fc=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0; H.val[950]=0; H.val[6947]=0
H.forward()
if H.val[37720]%9994531==0: H.val[8976]=H.val[37720]//9994531
H.forward()
cons=[287,1531,3081,7425,8273,8470,9708,13790,18030,27706,29926,30383,32138,37297,38740]
Cset=set(cons)
# nonzero atoms in consumers, with free inputs
seen=set()
for ai,poly in enumerate(atoms):
    if ateqs[ai]&Cset:
        val=ev(poly)
        if val!=0 and ai not in seen:
            seen.add(ai)
            vs=set()
            for m,c in poly: vs.update(m)
            fr=sorted(vs&H.freeinp)
            print('atom %d p?=%s eqs=%s free=%s :: %s'%(ai, val%p!=0, sorted(ateqs[ai]&Cset), fr, reprs[ai][:70]))
