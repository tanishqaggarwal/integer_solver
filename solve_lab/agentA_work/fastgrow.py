"""Equation-closure growth tracker.  At each level report atoms / equations / variables /
knobs / EXCLUDED (variables touching an atom outside the window) and whether the window is
still exactly affine (no monomial with >=2 knob factors).  The theorem becomes
UNCONDITIONAL when excluded = 0."""
import sys, collections, json, time; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
path=sys.argv[1]; LMAX=int(sys.argv[2]); EQBUD=int(sys.argv[3]) if len(sys.argv)>3 else 39033
v=L.load(path); fe=L.failing_eqs(L.all_atom_values(v))
print('%s failing=%d ; instance has %d atoms, %d equations, %d variables'%(
      path.split('/')[-1],len(fe),L.NA,L.NEQ,L.NVARS),flush=True)
A=set(a for e in fe for a in L.eq_atoms[e][2])
t0=time.time(); hist=[]
for lev in range(LMAX+1):
    R=set()
    for a in A: R|=set(L.atom2eq[a])
    V=set()
    for a in A: V|=L.avars[a]
    K=set(u for u in V if all(x in A for x in L.var_atoms[u]))
    excl=V-K
    nonlin=0
    for a in A:
        for m in L.polys[a]:
            c=sum(1 for u in m if u in K)
            if c>1: nonlin+=1; break
    print('L%-3d atoms=%-6d eqs=%-6d vars=%-6d knobs=%-6d EXCLUDED=%-6d nonlinear_atoms=%-5d [%.0fs]'%(
        lev,len(A),len(R),len(V),len(K),len(excl),nonlin,time.time()-t0),flush=True)
    hist.append({'L':lev,'atoms':len(A),'eqs':len(R),'vars':len(V),'knobs':len(K),
                 'excluded':len(excl),'nonlinear':nonlin})
    json.dump(hist,open('/home/user/integer_solver/solve_lab/agentA_work/growth.json','w'))
    if not excl:
        print('*** EXCLUDED SET IS EMPTY at L=%d -- the window is CLOSED and the bound is'%lev)
        print('    unconditional over all variables reachable from the residual. ***',flush=True)
        json.dump(sorted(A),open('/home/user/integer_solver/solve_lab/agentA_work/closed_window.json','w'))
        break
    if len(R)>EQBUD:
        print('equation budget %d exceeded; stopping'%EQBUD,flush=True); break
    A2=set()
    for e in R: A2|=set(L.eq_atoms[e][2])
    if A2==A:
        print('*** FIXED POINT at L=%d: the atom set is closed under equation closure.'%lev,flush=True)
        break
    A=A2
