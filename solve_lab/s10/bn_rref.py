"""bn_rref: sound strengthened peel using the ROW SPACE.

Every vector in the row space of M is a valid equality sum_a w_a t_a = 0.
With t >= 0, any such row whose coefficients are all one sign on the surviving
support forces every t in that support to 0.  Using RREF rows (and random
row-space combinations) gives a much stronger peel than the raw equations.
Iterate to a fixed point.
"""
import os, sys, json, collections, random, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

bools=B.bools_map()

def rows_for(S):
    idx={a:i for i,a in enumerate(S)}
    E=sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
    out=[]
    for e in E:
        m,sq,co=L.eq_atoms[e]
        r={}
        for a in co:
            if a in idx:
                w=co[a]*bools[a][1]           # coefficient on t_a
                if w: r[idx[a]]=Fraction(w)
        if r: out.append(r)
    return out

def rref(rows, n):
    R=[dict(r) for r in rows]; piv={}
    for ri in range(len(R)):
        r=R[ri]
        for c in sorted([c for c in list(r) if c in piv]):
            if r.get(c):
                f=r[c]/R[piv[c]][c]
                for k,v in R[piv[c]].items():
                    nv=r.get(k,Fraction(0))-f*v
                    if nv: r[k]=nv
                    else: r.pop(k,None)
        R[ri]=r
        if r: piv[min(r)]=ri
    pc=sorted(piv)
    for c in reversed(pc):
        r=R[piv[c]]; f=r[c]
        R[piv[c]]=r={k:v/f for k,v in r.items()}
        for c2 in pc:
            if c2>=c: continue
            r2=R[piv[c2]]
            if r2.get(c):
                g=r2[c]
                for k,v in r.items():
                    nv=r2.get(k,Fraction(0))-g*v
                    if nv: r2[k]=nv
                    else: r2.pop(k,None)
    return [R[piv[c]] for c in pc], pc

def peel(S, tries=4000, seed=0):
    S=sorted(S); rnd=random.Random(seed); rounds=0
    while True:
        rounds+=1
        n=len(S)
        rows=rows_for(S)
        rw,pc=rref(rows,n)
        kill=set()
        def check(r):
            if not r: return
            sg=set(1 if v>0 else -1 for v in r.values())
            if len(sg)==1: kill.update(r)
        for r in rw: check(r)
        # random row-space combinations
        if rw:
            for _ in range(tries):
                k=rnd.randint(2,min(4,len(rw)))
                cs=[rnd.choice(rw) for _ in range(k)]
                co=[Fraction(rnd.randint(-6,6)) for _ in range(k)]
                acc={}
                for c,r in zip(co,cs):
                    if not c: continue
                    for kk,v in r.items():
                        nv=acc.get(kk,Fraction(0))+c*v
                        if nv: acc[kk]=nv
                        else: acc.pop(kk,None)
                check(acc)
        # raw-equation peels too
        for r in rows:
            if len(r)==1: kill.update(r)
            else: check(r)
        if not kill:
            print(f'  fixed point after {rounds} rounds: {len(S)} atoms survive',flush=True)
            return S
        S=[S[i] for i in range(len(S)) if i not in kill]
        print(f'  round {rounds}: killed {len(kill)} -> {len(S)}',flush=True)
        if not S:
            print('  ALL ATOMS PEELED -> cone is TRIVIAL',flush=True)
            return []

if __name__=='__main__':
    S=json.load(open(os.path.join(HERE,'bn_keycomp.json')))['keycomp']
    print(f'row-space peel on the {len(S)}-atom maximal support',flush=True)
    t0=time.time()
    out=peel(S)
    print(f'RESULT: {len(out)} atoms survive the row-space peel ({time.time()-t0:.0f}s)',flush=True)
    json.dump({'survivors':out}, open(os.path.join(HERE,'bn_rref.json'),'w'))
