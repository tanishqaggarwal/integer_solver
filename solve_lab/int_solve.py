#!/usr/bin/env python3
"""Definitive integer solve of the 27 build_twist failing equations over the free-input handles,
using the wire (=p) product unpacking. Build integer Jacobian, solve A x = b over Z via HNF
(sympy). Apply, forward-eval, verify at equation level."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
import sympy
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
boolbits=set(json.load(open('boolbits.json'))['boolvars'])
override={22106:1, 16742:C2, 12186:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
handles=sorted(freeinp - boolbits - set(override)); hset=set(handles)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp or pinned[v]: ready[v]=True
gu=[0]*len(gates); using=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates):
    u=0
    for v in vids:
        if not ready[v]: u+=1
        using[v].append(gi)
    gu[gi]=u
definer={}; order=[]
q=deque(gi for gi in range(len(gates)) if gu[gi]==0)
while q:
    gi=q.popleft(); t,rhs,vids=gates[gi]
    if ready[t]: continue
    definer[t]=gi; order.append(t); ready[t]=True
    for gj in using[t]:
        gu[gj]-=1
        if gu[gj]==0: q.append(gj)
VAR=re.compile(r'x_(\d+)')
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
def inner_src(lhs):
    node=ast.parse(lhs,mode='eval').body
    while isinstance(node,ast.BinOp) and isinstance(node.op,ast.Mult):
        a,b=node.left,node.right
        ca=isinstance(a,ast.Constant) or (isinstance(a,ast.UnaryOp) and isinstance(a.operand,ast.Constant))
        cb=isinstance(b,ast.Constant) or (isinstance(b,ast.UnaryOp) and isinstance(b.operand,ast.Constant))
        if ca and not cb: node=b
        elif cb and not ca: node=a
        elif ast.unparse(a)==ast.unparse(b): node=a
        else: break
    return node
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
forward(); ns['v']=val
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"failing: {len(F)}", flush=True)
rootcode={i:compile(VAR.sub(r'v[\1]',ast.unparse(inner_src(lines[i].rsplit('=',1)[0]))),'<e>','eval') for i in F}
H=sorted(set().union(*[eqvars[i]&hset for i in F]))
base={i:eval(rootcode[i],ns) for i in F}
Jac=[[0]*len(H) for _ in F]
for j,h in enumerate(H):
    old=val[h]; val[h]=old+1; forward(); ns['v']=val
    for ri,i in enumerate(F): Jac[ri][j]=eval(rootcode[i],ns)-base[i]
    val[h]=old
forward()
# Solve A x = b over Z, A=Jac (m x n), b=-base. Use HNF: find unimodular U with A U = H (col HNF).
m=len(F); n=len(H)
Amat=Matrix(m,n,lambda i,j: Jac[i][j])
bvec=Matrix(m,1,lambda i,j: -base[F[i]])
# integer solve: stack [A | b]; a solution x exists iff b in column-lattice of A.
# Use the approach: compute HNF of A^T (rows=vars). Solve via the transform.
# Simpler robust: use sympy's linsolve for parametric Q-sol, then integer-search free params.
# --- Smith Normal Form with transforms: U A V = D ---
def snf_transform(Ain):
    m=len(Ain); n=len(Ain[0])
    Awrk=[row[:] for row in Ain]
    U=[[1 if i==j else 0 for j in range(m)] for i in range(m)]
    V=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    def swap_rows(i,j):
        Awrk[i],Awrk[j]=Awrk[j],Awrk[i]; U[i],U[j]=U[j],U[i]
    def swap_cols(i,j):
        for r in Awrk: r[i],r[j]=r[j],r[i]
        for r in V: r[i],r[j]=r[j],r[i]
    def addrow(i,j,f):  # row_i += f*row_j
        for k in range(n): Awrk[i][k]+=f*Awrk[j][k]
        for k in range(m): U[i][k]+=f*U[j][k]
    def addcol(i,j,f):  # col_i += f*col_j
        for r in range(m): Awrk[r][i]+=f*Awrk[r][j]
        for r in range(n): V[r][i]+=f*V[r][j]
    t=0
    for t in range(min(m,n)):
        # find a nonzero entry in submatrix >= (t,t) with min abs, bring to (t,t)
        while True:
            piv=None; best=None
            for i in range(t,m):
                for j in range(t,n):
                    if Awrk[i][j]!=0 and (best is None or abs(Awrk[i][j])<best): best=abs(Awrk[i][j]); piv=(i,j)
            if piv is None: return Awrk,U,V,t
            i0,j0=piv
            if i0!=t: swap_rows(i0,t)
            if j0!=t: swap_cols(j0,t)
            done=True
            for i in range(t+1,m):
                if Awrk[i][t]!=0:
                    f=-(Awrk[i][t]//Awrk[t][t]); addrow(i,t,f)
                    if Awrk[i][t]!=0: done=False
            for j in range(t+1,n):
                if Awrk[t][j]!=0:
                    f=-(Awrk[t][j]//Awrk[t][t]); addcol(j,t,f)
                    if Awrk[t][j]!=0: done=False
            if done:
                # ensure divisibility of remaining block
                bad=False
                for i in range(t+1,m):
                    for j in range(t+1,n):
                        if Awrk[i][j]%Awrk[t][t]!=0:
                            addrow(t,i,1); bad=True; break
                    if bad: break
                if not bad: break
    return Awrk,U,V,min(m,n)
D,U,V,rk=snf_transform([[Jac[i][j] for j in range(n)] for i in range(m)])
# c = U b
bb=[-base[F[i]] for i in range(m)]
c=[sum(U[i][k]*bb[k] for k in range(m)) for i in range(m)]
# solve D y = c
yv=[0]*n; feasible=True
for i in range(m):
    d=D[i][i] if i<n else 0
    if d==0:
        if c[i]!=0: feasible=False; break
    else:
        if c[i]%d!=0: feasible=False; break
        yv[i]=c[i]//d
# extra rows (i from n..m) require c[i]==0 already checked; free y (i>=rk) = 0
print(f"SNF integer solve: rank {rk}, feasible={feasible}", flush=True)
if not feasible:
    print("NO integer solution in the 48-handle subsystem -> diagnosing obstruction:")
    for i in range(m):
        d=D[i][i] if i<n else 0
        if d==0:
            if c[i]!=0: print(f"  row {i}: D=0 but c={c[i]} (2^{abs(c[i]).bit_length()}) -> rank/rhs mismatch")
        elif c[i]%d!=0:
            g=sympy.igcd(int(c[i]),int(d))
            print(f"  row {i}: D_ii=2^{abs(d).bit_length()} (p^{0}), c%D!=0, gcd(c,D)=2^{abs(g).bit_length()}, D/gcd factors p? {(abs(d)//abs(g))%p==0}")
    sys.exit(0)
# x = V y
xint=[sum(V[i][j]*yv[j] for j in range(n)) for i in range(n)]
for j,h in enumerate(H): val[h]=xint[j]
forward(); ns['v']=val
FF=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"INTEGER SOLVE (SNF) applied -> {len(lines)-len(FF)}/{len(lines)} ({len(FF)} fail): {FF[:20]}", flush=True)
if len(FF)<26:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('int_solved.json','w')); print("SAVED int_solved.json")
sys.exit(0)
solset=None  # (old parametric path disabled)
# find free symbols
freesyms=sorted(set().union(*[e.free_symbols for e in solset]), key=lambda s:s.name)
print(f"Q-solution parametric with {len(freesyms)} free integer params", flush=True)
# Each handle expression must be integer for integer free params. Collect denominators.
# Strategy: pick free params to clear denominators. Build congruences: for each handle expr,
# expr = (linear in free)/D must be integer => linear(free) ≡ num0 (mod D).
# Solve the congruence system with sympy (per-denominator).
from sympy import together, fraction, Rational, gcd as sgcd
# substitute free params = symbols; require each solset[i] in ZZ. Use solve_congruence chain.
# Build integer linear system for the fractional-part: for each i, expr_i has denominator d_i.
congr=[]  # (coeffs dict over freesyms, const, modulus)
for i,e in enumerate(solset):
    e=sympy.nsimplify(e); e=sympy.expand(e)
    num,den=fraction(together(e))
    den=int(den)
    if den==1: continue
    # num must be ≡0 mod den. num is linear in freesyms.
    poly=sympy.Poly(num, *freesyms) if freesyms else None
    row={}
    for s in freesyms:
        row[s]=int(poly.coeff_monomial(s)) if poly else 0
    const=int(poly.coeff_monomial(1)) if poly else int(num)
    congr.append((row,const,den))
print(f"{len(congr)} integrality congruences to satisfy", flush=True)
# Solve congruences greedily: treat as linear system mod lcm. Use sympy to solve.
# Set up: sum_s row[s]*s + const ≡ 0 (mod den). Solve for integer freesyms.
# Use sympy's linear_eq_to_matrix over each modulus is complex; instead brute force small: many
# free params -> set all but a few to 0, solve the rest via CRT/gcd per congruence.
freelist=list(freesyms)
# Report the denominators' p-power structure
denfactors=defaultdict(int)
for (_,_,den) in congr:
    d=den; k=0
    while d%p==0: d//=p; k+=1
    denfactors[(k,d)]+=1
print(f"congruence denominators (p-power k, residual d): {dict(denfactors)}", flush=True)
# Solve the mod-p reduction simultaneously over GF(p): sum_s row[s]*s ≡ -const (mod p) per congruence.
def invp(a): return pow(a%p, p-2, p)
grows=[]
for (row,const,den) in congr:
    grows.append(([row.get(s,0)%p for s in freelist], (-const)%p))
ncol=len(freelist)
Aug=[r[0][:]+[r[1]] for r in grows]
nrow=len(Aug); pr=0; pivc=[]
for c in range(ncol):
    sel=None
    for r in range(pr,nrow):
        if Aug[r][c]%p!=0: sel=r;break
    if sel is None: continue
    Aug[pr],Aug[sel]=Aug[sel],Aug[pr]
    iv=invp(Aug[pr][c]); Aug[pr]=[(x*iv)%p for x in Aug[pr]]
    for r in range(nrow):
        if r!=pr and Aug[r][c]%p!=0:
            f=Aug[r][c]; Aug[r]=[(Aug[r][k]-f*Aug[pr][k])%p for k in range(ncol+1)]
    pivc.append(c); pr+=1
    if pr>=nrow: break
incons_p=any(all(Aug[r][c]%p==0 for c in range(ncol)) and Aug[r][ncol]%p!=0 for r in range(nrow))
print(f"mod-p congruence system consistent: {not incons_p} (rank {len(pivc)}/{ncol})", flush=True)
assign={s:0 for s in freelist}
ok=not incons_p
if ok:
    # particular GF(p) solution (free non-pivot params = 0)
    solp={c:0 for c in range(ncol)}
    pivset=set(pivc)
    for r,c in enumerate(pivc): solp[c]=Aug[r][ncol]%p
    for c in range(ncol): assign[freelist[c]]=solp[c]
    print(f"mod-p solution found; free params set (nonzero: {sum(1 for c in range(ncol) if solp[c])})", flush=True)
if ok:
    # compute final handle values
    hv={}
    for i,e in enumerate(solset):
        v_=e.subs(assign)
        v_=sympy.nsimplify(v_)
        if v_.is_integer or (hasattr(v_,'q') and v_.q==1):
            hv[H[i]]=int(v_)
        else:
            print(f"  handle x_{H[i]} still non-integer: {v_}"); ok=False; break
    if ok:
        for h,x in hv.items(): val[h]=x
        forward(); ns['v']=val
        FF=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
        print(f"INTEGER SOLVE applied -> {len(lines)-len(FF)}/{len(lines)} ({len(FF)} fail): {FF[:20]}", flush=True)
        if len(FF)<26:
            json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('int_solved.json','w'))
            print("SAVED int_solved.json")
