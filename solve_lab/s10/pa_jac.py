"""Global sparse equation-Jacobian at the delivered witness, mod p."""
import os, sys, collections, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
P=2**256-2**32-977
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
vm=[x%P for x in v]
t0=time.time()
# atom jacobian rows
Jat={}
for a in range(L.NA):
    row={}
    for m,c in L.polys[a].items():
        for i,w in enumerate(m):
            t=c
            for j,z in enumerate(m):
                if j!=i: t=t*vm[z]%P
            row[w]=(row.get(w,0)+t)%P
    row={u:c for u,c in row.items() if c}
    if row: Jat[a]=row
print('atoms with nonzero grad:',len(Jat),'nnz',sum(len(r) for r in Jat.values()),f'{time.time()-t0:.0f}s')
# equation jacobian
K={}
for e,(mult,sq,co) in enumerate(L.eq_atoms):
    row=collections.defaultdict(int)
    for a,ca in co.items():
        r=Jat.get(a)
        if not r: continue
        for u,cu in r.items(): row[u]=(row[u]+ca*cu)%P
    row={u:c%P for u,c in row.items() if c%P}
    if row: K[e]=row
print('eqs with nonzero grad:',len(K),'nnz',sum(len(r) for r in K.values()),f'{time.time()-t0:.0f}s')
cols=set()
for r in K.values(): cols|=set(r)
print('live variables (columns):',len(cols))
av=L.all_atom_values(v)
r0={e: (L.eq_value(e,av))%P for e in range(L.NEQ)}
nzr=[e for e in r0 if r0[e]]
print('equations with nonzero residual mod p:',len(nzr),nzr)
json.dump({'K':{str(e):{str(u):str(c) for u,c in r.items()} for e,r in K.items()},
           'cols':sorted(cols),'nzr':nzr},open(os.path.join(HERE,'pa_K.json'),'w'))
print('saved',f'{time.time()-t0:.0f}s')
