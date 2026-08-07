import sys, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
P=2**256-2**32-977
chk=[a for a in range(L.NA) if a not in L.atom_out]
cls=collections.Counter(); ex=collections.defaultdict(list)
for a in chk:
    Pp=L.polys[a]; deg=max(len(m) for m in Pp)
    nv=len(L.avars[a]); nt=len(Pp)
    # boolean: x^2 - x
    key=None
    if deg==2 and nv==1:
        u=next(iter(L.avars[a]))
        if Pp.get((u,u))==1 and Pp.get((u,))==-1 and nt==2: key='BOOL'
    if key is None and deg==1 and nv==1:
        key='PIN'  # c*x + k
    if key is None and deg==4: key='SQUARE4'
    if key is None: key=f'deg{deg}_nv{nv if nv<4 else "many"}_nt{nt if nt<5 else "many"}'
    cls[key]+=1
    if len(ex[key])<3: ex[key].append((a,L.atom_src[a][:150]))
for k,c in cls.most_common(25):
    print(f'{c:6d}  {k}')
print()
for k in ['PIN','BOOL','SQUARE4']:
    for a,s in ex[k]: print(k,a,s)
