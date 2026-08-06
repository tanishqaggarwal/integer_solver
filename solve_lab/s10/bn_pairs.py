"""bn_pairs: EXACT minimum equation-cost of a two-atom boolean carrier.

For a pair (a,b) with values t_a,t_b>0, a shared equation e cancels iff
   coeff(e,a)*c_a*t_a + coeff(e,b)*c_b*t_b = 0
i.e. t_a/t_b = -coeff(e,b)*c_b / (coeff(e,a)*c_a) =: r_e   (must be > 0).
All shared equations with the SAME positive ratio cancel simultaneously.
So   mincost(a,b) = |E(a) u E(b)| - max_r #{e shared : r_e = r > 0}.
Then apply the k(k-1) realisability filter to that ratio.

Chunked + checkpointed to JSONL.  usage: bn_pairs.py START END
"""
import os, sys, json, collections, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'
sys.path.insert(0, HERE)
import bn_lib as B
import lib as L, tools as T, ad

OUT=os.path.join(HERE,'bn_pairs.jsonl')
bools=B.bools_map()
FREEB=sorted(a for a in bools if bools[a][0] in B.FREESET)
FSET=set(FREEB)
downstream={}                    # measured downstream damage per atom
for l in open(os.path.join(HERE,'bn_sweep.jsonl')):
    r=json.loads(l); downstream[r['atom']]=39026-r['score']-len(L.atom2eq[r['atom']])

# pairs sharing at least one equation
e2a=collections.defaultdict(list)
for a in FREEB:
    for e in L.atom2eq[a]: e2a[e].append(a)
pairs=set()
for e,As in e2a.items():
    if len(As)>200: continue
    for i in range(len(As)):
        for j in range(i+1,len(As)): pairs.add((As[i],As[j]))
pairs=sorted(pairs)
print('free boolean atom pairs sharing >=1 equation:',len(pairs),flush=True)

TRI={k*(k-1):k for k in range(1,4000)}
TRIS=sorted(t for t in TRI if t>0)
def realisable(r):
    """is there t_a=k(k-1), t_b=m(m-1) with t_a/t_b == r?"""
    p,q=r.numerator,r.denominator
    for tb in TRIS[:400]:
        ta=Fraction(tb*p,q)
        if ta.denominator==1 and int(ta) in TRI:
            return (int(ta),tb)
    return None

S=int(sys.argv[1]); E=min(int(sys.argv[2]),len(pairs))
done=set()
if os.path.exists(OUT):
    for l in open(OUT):
        try: done.add(json.loads(l)['k'])
        except Exception: pass
f=open(OUT,'a'); t0=time.time(); best=(10**9,)
for i in range(S,E):
    a,b=pairs[i]; k=f'{a},{b}'
    if k in done: continue
    Ea=set(L.atom2eq[a]); Eb=set(L.atom2eq[b])
    sh=Ea&Eb
    ca,cb=bools[a][1],bools[b][1]
    cnt=collections.Counter()
    for e in sh:
        m,sq,co=L.eq_atoms[e]
        num=-co[b]*cb; den=co[a]*ca
        r=Fraction(num,den)
        if r>0: cnt[r]+=1
    tot=len(Ea|Eb)
    if cnt:
        r,n=cnt.most_common(1)[0]
        cost=tot-n
        rl=realisable(r)
    else:
        r=None; n=0; cost=tot; rl=None
    dsn=downstream.get(a,0)+downstream.get(b,0)
    rec={'k':k,'a':a,'b':b,'cost':cost,'shared':len(sh),'cancel':n,
         'ratio':(str(r) if r else None),'real':rl,'down':dsn,'tot_est':cost+dsn}
    f.write(json.dumps(rec)+'\n')
    if rl and cost+dsn<best[0]:
        best=(cost+dsn,rec); print('  cand',rec,flush=True)
f.flush(); f.close()
print(f'[{S},{E}) done {time.time()-t0:.0f}s best={best[0]}',flush=True)
