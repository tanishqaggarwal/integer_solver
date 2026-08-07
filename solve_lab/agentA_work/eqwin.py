"""Equation-level windows around the deliverable.
A_0 = atoms of the failing equations.  Equation-closure step: R = eqs(A), A := atoms(R).
This explicitly ADMITS the foreign atoms of every modelled equation as free-to-be-nonzero,
i.e. it opens exactly the cancellation freedom that the atom-closure windows suppressed.
Report, at each level: atoms, equations, unfiltered knobs, knobs after the linearity
filter, how many atoms are nonlinear in the unfiltered knobs, and how many of those are
perfect squares (rescuable to an affine row)."""
import sys, collections, math, json; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import pick_knobs, build, qsolve
P=env.P
path=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
LEV=int(sys.argv[2]) if len(sys.argv)>2 else 3
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av)
print('%s score=%d failing=%d'%(path.split('/')[-1],L.NEQ-len(fe),len(fe)),flush=True)
def knobpoly(a,K,v):
    ki={u:i for i,u in enumerate(K)}
    out=collections.defaultdict(int)
    for m,c in L.polys[a].items():
        t=c; mono=[]
        for u in m:
            if u in ki: mono.append(ki[u])
            else: t*=v[u]
        if t: out[tuple(sorted(mono))]+=t
    return {k:c for k,c in out.items() if c}
def is_square(kp):
    """is the knob-polynomial a perfect square of an affine form?"""
    vs=sorted(set(i for m in kp for i in m))
    if any(len(m)>2 for m in kp): return False
    q={}
    for i in vs:
        c=kp.get((i,i),0)
        if c<0: return False
        s=math.isqrt(c)
        if s*s!=c: return False
        q[i]=s
    c0=kp.get((),0)
    if c0<0: return False
    s0=math.isqrt(c0)
    if s0*s0!=c0: return False
    if not vs: return False
    ref=vs[0]
    sg={ref:1}
    for i in vs[1:]:
        cij=kp.get(tuple(sorted((ref,i))),0)
        sg[i]=1 if cij>0 else -1
    for sc0 in (1,-1):
        Q={i:sg[i]*q[i] for i in vs}; Qc=sc0*s0
        chk=collections.defaultdict(int); chk[()]=Qc*Qc
        for i in vs:
            chk[(i,)]+=2*Qc*Q[i]
            for j in vs: chk[tuple(sorted((i,j)))]+=Q[i]*Q[j]
        if all(chk.get(m,0)==kp.get(m,0) for m in set(chk)|set(kp)): return True
    return False
A=set(a for e in fe for a in L.eq_atoms[e][2])
for lev in range(LEV+1):
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    Kraw=sorted(u for u in set(u for a in A for u in L.avars[a]) if all(x in A for x in L.var_atoms[u]))
    nl=[]; sq=0
    for a in A:
        kp=knobpoly(a,Kraw,v)
        if any(len(m)>1 for m in kp):
            nl.append(a)
            if is_square(kp): sq+=1
    Kf=pick_knobs(v,A)
    print('lev%d  atoms=%-5d eqs=%-5d knobs_raw=%-4d knobs_linear=%-4d nonlinear_atoms=%-4d of which perfect squares=%d'%(
        lev,len(A),len(R),len(Kraw),len(Kf),len(nl),sq),flush=True)
    if lev<LEV:
        A=set(a for e in R for a in L.eq_atoms[e][2])   # EQUATION closure: admit foreign atoms
