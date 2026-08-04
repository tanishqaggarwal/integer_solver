import os,sys
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import heal_harness as H, json
p=H.p
atoms=[]; reprs={}
with open('atoms/poly_atoms.jsonl') as f:
    for i,line in enumerate(f):
        d=json.loads(line); atoms.append([(tuple(m),c) for m,c in d['poly']]); reprs[i]=d.get('repr','')
def atomnz(v):
    nz=[]
    for ai,poly in enumerate(atoms):
        s=0
        for m,c in poly:
            t=c
            for var in m: t*=v[var]
            s+=t
        if s!=0: nz.append((ai,s))
    return nz
def setup(overrides,absorb=True):
    vA=H.loadd('best_agentA_39022.json')
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    for k,val in overrides.items(): H.val[k]=val
    H.forward()
    return vA
if __name__=='__main__':
    setup({4287:1})
    t4432=H.val[4432]; t7068=H.val[7068]
    H.val[8731]=t4432; H.val[9118]=t7068; H.val[9413]=0; H.val[17325]=0
    H.forward()
    F=H.fails(); print('fails=',len(F))
    nz=atomnz(H.val)
    print('nonzero atoms=',len(nz))
    for ai,s in nz:
        print(f'atom {ai}: zeroModP={s%p==0} val%p={s%p if s%p!=0 else 0} :: {reprs[ai][:85]}')
