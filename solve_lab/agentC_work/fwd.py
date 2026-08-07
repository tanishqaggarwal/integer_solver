import sys, collections, pickle, json, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentC_work/supp3.pkl','rb'))
outs=D['outs']; free=D['free']; topo=D['topo']; csupp=D['csupp']
P=2**256-2**32-977
# precompute for each defined var: (coeff, [(coef, monomial-without-out)]) for linear-in-out solve
lin={}
bad=0
for v,a in outs.items():
    c=0; rest=[]
    ok=True
    for m,cc in L.polys[a].items():
        k=m.count(v)
        if k>1: ok=False; break
        if k==0: rest.append((cc,m))
        else: 
            mm=tuple(w for w in m if w!=v)
            if mm: ok=False; break   # out multiplied by another var
            c+=cc
    if not ok or c==0: bad+=1; continue
    lin[v]=(c,rest)
print('gates solvable linearly for output:',len(lin),'bad:',bad)
cs=collections.Counter(abs(c) for c,_ in lin.values())
print('|coeff| of output:',cs.most_common(5))
def forward(vals):
    v=vals
    ndiv=0
    for u in topo:
        if u in lin:
            c,rest=lin[u]
            s=0
            for cc,m in rest:
                t=cc
                for w in m: t*=v[w]
                s+=t
            if s % c: ndiv+=1; v[u]=0
            else: v[u]=-s//c
    return ndiv
v=[0]*L.NVARS
nd=forward(v)
print('nondivisible gates at all-zero:',nd)
av=L.all_atom_values(v)
chk=[a for a in range(L.NA) if a not in L.atom_out]
nzchk=[a for a in chk if av[a]!=0]
print('nonzero checks at all-zero-free:',len(nzchk))
zs=[a for a in chk if not csupp[a]]
print('constant checks nonzero:',sum(1 for a in zs if av[a]!=0))
fails=L.failing_eqs(av)
print('score at all-zero-free:',L.NEQ-len(fails))
# also nonzero gate atoms
nzg=[a for a in range(L.NA) if a in L.atom_out and av[a]!=0]
print('nonzero gate atoms:',len(nzg))
json.dump({'nzchk':nzchk[:50]},open('/dev/null','w'))
pickle.dump({'lin':lin},open('/home/user/integer_solver/solve_lab/agentC_work/lin.pkl','wb'))
