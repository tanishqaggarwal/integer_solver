#!/usr/bin/env python3
"""Quadratic-core solve with CONFLICT-FREE null space: iteratively freeze handles whose product-
partners also move (2nd-order breakers), recompute achievable (da,dc), re-solve S=T=0, until the
applied move breaks no wiring eq. Then set quotient handles and verify; save if >39013."""
import json, time, ast, re, sys
from agentC_common import (p, gates, order, definer, gcode, forward, partial_forward, downstream_ks,
                           val, freeinp, ns, lines, eqcode, eqvars, load_best, CORE, posof, NVARS,
                           pinned, rootcode_of, inv, is_qr, sqrt_mod, C1, C2)
import agentC_poly as Ply
from collections import defaultdict

ROOTSEL = int(sys.argv[1]) if len(sys.argv)>1 else 1   # 0=regime1, 1=regime2
best=load_best(); forward()
gate_defs={t:(rhs,vids) for t,rhs,vids in gates}
_fcmemo={}
def freecone(root):
    if root in _fcmemo: return _fcmemo[root]
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    _fcmemo[root]=leaves; return leaves
DEEP=[3558,29322,33469,27713,1326]
deepcone=set()
for d in DEEP+[35389,6671]: deepcone|=freecone(d)
a_id,b_id,c_id=12186,14853,16742
CONTROLS=[a_id,b_id,c_id]
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
def affected_eqs(h):
    aff=set(eqbyvar.get(h,()))
    for k in downstream_ks(h): aff|=eqbyvar.get(order[k],set())
    return aff
# closure
H=set(CONTROLS); cons=set()
for _hop in range(8):
    newcons=set()
    for h in H: newcons|=(affected_eqs(h)-F0)
    cons=newcons
    comp=set()
    for i in cons: comp|=(eqvars[i]&freeinp)
    comp-=deepcone; comp|=set(CONTROLS)
    if comp<=H: break
    H=comp
H=sorted(H); Hidx={h:i for i,h in enumerate(H)}; NH=len(H); cons=sorted(cons)
print(f"handles={NH}, constraints={len(cons)}", flush=True)

# product-partner map among handles (within cons) for conflict detection
VAR=re.compile(r'x_(\d+)')
def product_pairs(expr):
    node=ast.parse(expr,mode='eval').body; pairs=[]
    def vids(n): return set(int(m) for m in VAR.findall(ast.unparse(n)))
    def walk(n):
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Mult):
                la=vids(n.left); lb=vids(n.right)
                if la and lb: pairs.append((la,lb))
            walk(n.left); walk(n.right)
        elif isinstance(n,ast.UnaryOp): walk(n.operand)
    walk(node); return pairs
Hset=set(H)
partner=defaultdict(set)
prod_eqs=[]   # (eqidx, [(coneA_handles, coneB_handles)])
for i in cons:
    pl=[]
    for (la,lb) in product_pairs(lines[i].rsplit('=',1)[0]):
        ca=set(); cb=set()
        for v in la: ca|=freecone(v)
        for v in lb: cb|=freecone(v)
        cah=ca&Hset; cbh=cb&Hset
        if cah and cbh:
            pl.append((cah,cbh))
            for x in cah: partner[x]|=cbh-{x}
            for x in cbh: partner[x]|=cah-{x}
    if pl: prod_eqs.append((i,pl))
print(f"constraint eqs with handle-products: {len(prod_eqs)}", flush=True)

# Jacobian
base_root={i:eval(rootcode_of(i),ns)%p for i in cons}
dks={h:downstream_ks(h) for h in H}
consset=set(cons)
aff_h={}
for h in H:
    a=set(eqbyvar.get(h,()))&consset
    for k in dks[h]: a|=(eqbyvar.get(order[k],())&consset)
    aff_h[h]=sorted(a)
Jcols=[]
for h in H:
    o=val[h]; val[h]=o+1; partial_forward(dks[h])
    col={}
    for i in aff_h[h]:
        d=(eval(rootcode_of(i),ns)-base_root[i])%p
        if d: col[i]=d
    Jcols.append(col); val[h]=o; partial_forward(dks[h])
rows_used=sorted(set().union(*[set(c) for c in Jcols])) if Jcols else []
def build_null(extra_zero):
    """null space of J plus rows forcing handle in extra_zero to 0. Returns list of null vecs."""
    ridx={r:i for i,r in enumerate(rows_used)}; M=len(rows_used)
    Jm=[[0]*NH for _ in range(M)]
    for hj,col in enumerate(Jcols):
        for r,vv in col.items(): Jm[ridx[r]][hj]=vv%p
    for h in extra_zero:
        row=[0]*NH; row[Hidx[h]]=1; Jm.append(row)
    # rref
    Mm=[row[:] for row in Jm]; m=len(Mm); r=0; piv=[]
    for c in range(NH):
        sel=None
        for i in range(r,m):
            if Mm[i][c]%p: sel=i;break
        if sel is None: continue
        Mm[r],Mm[sel]=Mm[sel],Mm[r]; ivv=inv(Mm[r][c]); Mm[r]=[(x*ivv)%p for x in Mm[r]]
        for i in range(m):
            if i!=r and Mm[i][c]%p:
                f=Mm[i][c]; Mm[i]=[(Mm[i][k]-f*Mm[r][k])%p for k in range(NH)]
        piv.append(c); r+=1
        if r>=m: break
    pivset=set(piv); free_cols=[c for c in range(NH) if c not in pivset]
    Null=[]
    for f in free_cols:
        n=[0]*NH; n[f]=1
        for ri,c in enumerate(piv): n[c]=(-Mm[ri][f])%p
        Null.append(n)
    return Null,len(piv)
def rref_small(Mrows,ncol):
    Mm=[row[:] for row in Mrows]; m=len(Mm); r=0; piv=[]
    for c in range(ncol):
        sel=None
        for i in range(r,m):
            if Mm[i][c]%p: sel=i;break
        if sel is None: continue
        Mm[r],Mm[sel]=Mm[sel],Mm[r]; ivv=inv(Mm[r][c]); Mm[r]=[(x*ivv)%p for x in Mm[r]]
        for i in range(m):
            if i!=r and Mm[i][c]%p:
                f=Mm[i][c]; Mm[i]=[(Mm[i][k]-f*Mm[r][k])%p for k in range(ncol)]
        piv.append(c); r+=1
        if r>=m: break
    return Mm[:r],piv

ia,ib,ic=Hidx[a_id],Hidx[b_id],Hidx[c_id]
x29322=val[29322]%p; x1326=val[1326]%p; x27713=val[27713]%p; x33469=val[33469]%p; x3558=val[3558]%p
T0=val[6671]%p
def solve_quadratic(projbasis, nullbasis):
    """projbasis: list of [da,db,dc]; nullbasis: full moves. Returns list of (da,dc,fullmove_dict)."""
    # need db component 0 in achievable? we solve generally: param t in R^dimV.
    dimV=len(projbasis)
    # Express (da,db,dc) = sum t_j projbasis[j]; require the resulting deep vars solve S=T=0.
    # We only have da,dc affecting S,T (db enters x29322=b-a,x33469=a+b+K33). Keep general:
    # da_t, db_t, dc_t linear in t. Substitute into S,T -> polynomial system in t (dimV vars).
    # For our case projbasis columns; handle dimV up to 3 by picking regime.
    return dimV
# Iteration: start with no frozen, detect conflicts from actual applied move, freeze, repeat.
frozen=set()
def moved_from(nullbasis, projbasis, da, dc):
    # map (da,dc) to full move using basis where projbasis rows are [1,0,0]/[0,0,1] style
    pass

# --- main loop ---
K33=(C1 + 97553848499418123410591666447050222001188385549510401465815187079080512838891)%p
def get_basis(frozen):
    Null,rank=build_null(frozen)
    Proj=[[n[ia]%p,n[ib]%p,n[ic]%p] for n in Null]
    _,pp=rref_small(Proj,3); dimV=len(pp)
    # choose basis vectors with independent projections
    Ptmp=[]; idx=[]
    for j,pj in enumerate(Proj):
        if len(rref_small(Ptmp+[pj],3)[1])>len(rref_small(Ptmp,3)[1]):
            Ptmp.append(pj); idx.append(j)
            if len(idx)==dimV: break
    return Null,[Proj[j] for j in idx],[Null[j] for j in idx],dimV
def solve_and_move(projb, nullb):
    """Solve S=T=0 given projection basis (want to control da,dc). Return list of (fullmove, info)."""
    dimV=len(projb)
    # Build linear maps: da(t),db(t),dc(t)
    def comp(coord): return [projb[j][coord] for j in range(dimV)]
    Da,Db,Dc=comp(0),comp(1),comp(2)
    # deep residues as functions of t (linear): need x29322(t)=x29322 + db - da ; etc.
    # We'll pick regime via ROOTSEL by parameterizing. Simplify: if we can get da,dc independent and
    # db=0 (or expressible), reduce to (da,dc). General handling: solve S=T=0 by treating t as unknown
    # with elimination when dimV==2.
    sols=[]
    if dimV==2:
        # t=(t0,t1). da=Da0 t0+Da1 t1, etc.
        # deep vars linear in t:
        # x29322_t = x29322 + (Db-Da).t ; x1326_t = x1326 + Da.t ; x33469_t = x33469 + (Da+Db).t
        # x27713_t = x27713 + Dc.t ; x3558_t = x3558 - Dc.t
        # unknowns t0,t1. Solve T=0 (deg2) & S=0(deg3). Use resultant via Ply on t0 after solving?
        # Represent everything as polynomials in t0 with coeff polynomials in t1? Do double elimination.
        # Simpler: sample-free direct — reduce using T then S via 2D Groebner-ish elimination by
        # treating as polynomials in t1 with coefficients in F_p[t0]; compute resultant.
        import agentC_bipoly as BP
        # coefficients (constant, t0, t1)
        def aff(base, v):  # base + v0*t0 + v1*t1  -> BP poly
            return BP.mk(base, v[0], v[1])
        dA=aff(0,Da); dB=aff(0,Db); dC=aff(0,Dc)
        X29=BP.add(BP.const(x29322), BP.sub(dB,dA))
        X13=BP.add(BP.const(x1326), dA)
        X33=BP.add(BP.const(x33469), BP.add(dA,dB))
        X27=BP.add(BP.const(x27713), dC)
        X35=BP.sub(BP.const(x3558), dC)
        S=BP.sub(BP.mul(X33,BP.mul(X29,X29)), BP.mul(X35,X35))
        T=BP.sub(BP.mul(X27,X29), BP.mul(X35,X13))
        for (t0v,t1v) in BP.solve2(S,T):
            tv=[t0v,t1v]
            da=sum(Da[j]*tv[j] for j in range(2))%p
            db=sum(Db[j]*tv[j] for j in range(2))%p
            dc=sum(Dc[j]*tv[j] for j in range(2))%p
            fm={h:(tv[0]*nullb[0][Hidx[h]]+tv[1]*nullb[1][Hidx[h]])%p for h in H}
            sols.append((da,db,dc,fm))
    return sols

for it in range(12):
    Null,projb,nullb,dimV=get_basis(frozen)
    print(f"iter {it}: frozen={len(frozen)}, nullity={len(Null)}, achievable dim={dimV}, projbasis={projb}", flush=True)
    if dimV<2:
        print("  achievable dim < 2, cannot solve conic; stop"); break
    sols=solve_and_move(projb,nullb)
    print(f"  quadratic solutions: {len(sols)}", flush=True)
    applied=False
    for (da,db,dc,fm) in sols:
        moved=set(h for h in H if fm[h]!=0)
        # detect conflicts: product eq with both cones having moved handles
        confl=set()
        for (i,pl) in prod_eqs:
            for (cah,cbh) in pl:
                if (cah&moved) and (cbh&moved):
                    confl|=((cah|cbh)&moved)
        # allow controls to be 'moved' (their products are core) -> exclude control-only products
        confl-=set(CONTROLS)
        if confl:
            continue  # this solution has conflicts; try next / will freeze below
        # apply and measure
        snap={h:val[h] for h in H}
        for h in H: val[h]=val[h]+fm[h]
        forward(); ns['v']=val
        if val[11150]%p==0: val[30317]=-(val[11150]//p)
        if (537773*val[37758])%p==0: val[2936]=(537773*val[37758])//p
        if val[25739]%(6672769*p)==0: val[5146]=val[25739]//(6672769*p)
        forward(); ns['v']=val
        F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
        core=[i for i in F if i in CORE]; nc=[i for i in F if i not in CORE]
        print(f"  APPLIED da/dc: sat={len(lines)-len(F)}/{len(lines)} core={len(core)} noncore={len(nc)} S%p={val[35389]%p==0} T%p={val[6671]%p==0}")
        if len(nc)==0:
            print(f"  *** all non-core satisfied! core={sorted(core)}")
            json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open(f'agentC_state_{len(lines)-len(F)}.json','w'))
            print(f"  saved agentC_state_{len(lines)-len(F)}.json")
        # keep applied for reporting; if worse restore
        if len(F) < 20:
            print("  IMPROVEMENT retained"); applied=True; break
        for h in H: val[h]=snap[h]
        forward(); ns['v']=val
        nc_meas=nc
    if applied: break
    # no conflict-free solution: find conflicts from FIRST solution and freeze
    if not sols: print("  no solutions"); break
    da,db,dc,fm=sols[0]
    moved=set(h for h in H if fm[h]!=0)
    confl=set()
    for (i,pl) in prod_eqs:
        for (cah,cbh) in pl:
            if (cah&moved) and (cbh&moved):
                confl|=((cah|cbh)&moved)
    confl-=set(CONTROLS)
    newf=confl-frozen
    if not newf:
        print("  no new conflicts to freeze but solution still breaks; measuring & stopping")
        break
    print(f"  freezing {len(newf)} conflict handles: {sorted(newf)[:20]}")
    frozen|=newf
