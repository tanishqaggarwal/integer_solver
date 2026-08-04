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
# for each consumer eq, list ALL its atoms: value, and free slacks (free var with unit-ish coeff, non-p)
for eq in [287,8273,13790,18030]:
    print('=== eq %d ==='%eq)
    for ai,poly in enumerate(atoms):
        if eq in ateqs[ai]:
            val=ev(poly)
            vs=set()
            for m,c in poly: vs.update(m)
            fr=sorted(vs&H.freeinp)
            # identify unit-coeff free slacks (free var appearing linearly with small coeff)
            slacks=[]
            for var in fr:
                cc=0; nl=False
                for m,c in poly:
                    if m.count(var)==1:
                        t=c
                        for v2 in m:
                            if v2!=var: t*=H.val[v2]
                        cc+=t
                    elif m.count(var)>1: nl=True
                if not nl and cc!=0:
                    gran='p' if cc%p==0 else ('unit' if abs(cc)==1 else str(cc%p)[:6])
                    slacks.append((var,gran))
            tag='NZ' if val!=0 else '..'
            if val!=0 or slacks:
                print('  [%s] atom %d free=%s slacks=%s'%(tag,ai,fr[:6],slacks[:6]))
