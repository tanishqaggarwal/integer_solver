import sys, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentC_work/supp3.pkl','rb'))
LN=pickle.load(open('/home/user/integer_solver/solve_lab/agentC_work/lin.pkl','rb'))['lin']
csupp=D['csupp']; topo=D['topo']
P=2**256-2**32-977
def forward(v):
    for u in topo:
        if u in LN:
            c,rest=LN[u]; s=0
            for cc,m in rest:
                t=cc
                for w in m: t*=v[w]
                s+=t
            v[u]=-s//c if s%c==0 else 0
    return v
v=forward([0]*L.NVARS)
av=L.all_atom_values(v)
chk=[a for a in range(L.NA) if a not in L.atom_out]
nz=[a for a in chk if av[a]!=0]
for a in nz:
    print('== a%d  val=%d  (val mod p = %d)'%(a,av[a],av[a]%P))
    print('   free-supp size',len(csupp[a]),'  eqs',len(L.atom2eq.get(a,{})))
    print('   src:',L.atom_src[a][:400])
    print('   supp:',sorted(csupp[a])[:40])
fails=L.failing_eqs(av)
print('failing eqs',len(fails),fails[:30])
