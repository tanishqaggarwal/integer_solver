"""bn_kernel: exact rank / nullspace of boolean-atom incidence blocks, plus the
k(k-1) realisability filter.

Atom value of boolean atom a on var u is  val_a = c_a * (x_u^2 - x_u) = c_a * t
with t = k(k-1) in {0,2,6,12,20,30,...}  (ALWAYS >= 0, always even,
t/2 triangular).  So writing t_a >= 0 the equation constraint becomes
   sum_a coeff(e,a) * c_a * t_a = 0     for every e in E(S),  t >= 0.

SIGN PEEL: if within an equation every surviving atom has the same sign of
coeff*c_a, all those t_a are forced to 0.  Iterate to a fixed point.
"""
import os, sys, json, collections
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools = B.bools_map()
BA=set(bools)

def sign_peel(cand, verbose=True):
    cand=set(cand); it=0
    while True:
        it+=1
        eqs=collections.defaultdict(list)
        for a in cand:
            for e,co in L.atom2eq[a].items(): eqs[e].append((a, co*bools[a][1]))
        kill=set()
        for e,As in eqs.items():
            sgn=set((1 if c>0 else -1) for _,c in As)
            if len(sgn)==1:                 # all same sign, t>=0 => all zero
                kill.update(a for a,_ in As)
        if not kill: break
        cand-=kill
        if verbose: print(f'   sign-peel round {it}: killed {len(kill)} -> {len(cand)}')
    return cand

def build(S):
    S=sorted(S)
    idx={a:i for i,a in enumerate(S)}
    E=sorted(set().union(*[set(L.atom2eq[a]) for a in S])) if S else []
    rows=[]
    for e in E:
        m,sq,co=L.eq_atoms[e]
        r={idx[a]: co[a] for a in co if a in idx}
        if r: rows.append(r)
    return S, E, rows

def nullspace(S, rows):
    """exact rational nullspace of the sparse row system.  Returns (rank, basis)"""
    n=len(S)
    R=[dict((k,Fraction(v)) for k,v in r.items()) for r in rows]
    piv={}          # col -> row index
    order=[]
    for ri,r in enumerate(R):
        # reduce against existing pivots
        for c in sorted([c for c in r if c in piv]):
            if c not in r or r[c]==0: continue
            f=r[c]/R[piv[c]][c]
            for k,vv in R[piv[c]].items():
                nv=r.get(k,Fraction(0))-f*vv
                if nv: r[k]=nv
                else: r.pop(k,None)
        r={k:v for k,v in r.items() if v}
        R[ri]=r
        if r:
            c=min(r); piv[c]=ri; order.append(c)
    rank=len(piv)
    # back-substitute to reduced form
    pc=sorted(piv)
    for c in reversed(pc):
        ri=piv[c]; r=R[ri]
        f=r[c]
        R[ri]=r={k:v/f for k,v in r.items()}
        for c2 in pc:
            if c2>=c: continue
            r2=R[piv[c2]]
            if c in r2 and r2[c]:
                g=r2[c]
                for k,vv in r.items():
                    nv=r2.get(k,Fraction(0))-g*vv
                    if nv: r2[k]=nv
                    else: r2.pop(k,None)
    freec=[c for c in range(n) if c not in piv]
    basis=[]
    for fc in freec:
        vec=[Fraction(0)]*n; vec[fc]=Fraction(1)
        for c in pc:
            r=R[piv[c]]
            if fc in r: vec[c]=-r[fc]
        basis.append(vec)
    return rank, basis, freec

def istri(t):
    """is t == k(k-1) for some integer k?  (t>=0, even, 1+4t a perfect square)"""
    if t<0 or t%2: return False
    m=1+4*t
    r=int(m**0.5)
    while r*r>m: r-=1
    while (r+1)*(r+1)<=m: r+=1
    return r*r==m

def report(name, S):
    S,E,rows = build(S)
    print(f'\n=== {name}: |S|={len(S)} |E(S)|={len(E)} defic={len(E)-len(S)} ===')
    if not S: return
    rank, basis, freec = nullspace(S, rows)
    print(f'   exact rank over Q = {rank}, nullspace dim = {len(S)-rank}')
    if not basis:
        print('   NO KERNEL: every boolean atom in this block is forced to zero.')
        return
    # sign structure of the kernel
    print(f'   kernel basis vectors: {len(basis)}')
    for bi,vec in enumerate(basis[:3]):
        supp=[(S[i], vec[i]) for i in range(len(S)) if vec[i]]
        pos=sum(1 for a,x in supp if x*bools[a][1]>0)
        neg=len(supp)-pos
        print(f'     v{bi}: support {len(supp)}  t>0:{pos} t<0:{neg}')
    return S,E,rows,rank,basis

if __name__=='__main__':
    core=set(json.load(open(os.path.join(HERE,'bn_core.json'))))
    fcore=set(json.load(open(os.path.join(HERE,'bn_fcore.json'))))
    d=json.load(open(os.path.join(HERE,'bn_defic.json')))
    blk=set(d['all']['S'])

    print('### SIGN PEEL (t = x(x-1) >= 0 always) ###')
    print(' -- all boolean atoms --')
    sp_all = sign_peel(BA)
    print(f'   survivors {len(sp_all)}')
    print(' -- zero-peel core (2172) --')
    sp_core = sign_peel(core)
    print(f'   survivors {len(sp_core)}')
    print(' -- free-var core (241) --')
    sp_f = sign_peel(fcore)
    print(f'   survivors {len(sp_f)}')
    print(' -- deficiency -29 block (376) --')
    sp_b = sign_peel(blk)
    print(f'   survivors {len(sp_b)}')

    report('free-var zero-core (241)', fcore)
    report('deficiency -29 block (376)', blk)
    json.dump({'sp_all':sorted(sp_all),'sp_core':sorted(sp_core),
               'sp_free':sorted(sp_f),'sp_blk':sorted(sp_b)},
              open(os.path.join(HERE,'bn_sign.json'),'w'))
