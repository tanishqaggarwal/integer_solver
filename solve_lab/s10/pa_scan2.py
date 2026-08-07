"""Deficiency scan over PRIMITIVE atoms only (bundles excluded: a bundle's value is
determined by the primitives so it adds +1 to c, cancelling its -1 combinatorial gain).
Seed = each primitive atom; then greedy equation growth inside a local pool."""
import os, sys, collections, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
EQS=[frozenset(L.atom2eq[a]) for a in range(L.NA)]
PRIM=[a for a in range(L.NA) if len(EQS[a])>1]
isprim=[len(EQS[a])>1 for a in range(L.NA)]
eqat=[[a for a in co if isprim[a]] for m,sq,co in L.eq_atoms]
t0=time.time()
res=[]
for a in PRIM:
    Ea=EQS[a]
    cnt=collections.Counter()
    for e in Ea:
        for b in eqat[e]: cnt[b]+=1
    ins=[b for b,c in cnt.items() if c==len(EQS[b])]
    res.append((len(ins)-len(Ea),a,len(Ea),tuple(sorted(ins))))
res.sort(key=lambda t:-t[0])
print('TOP seeds, f = |A(eqs(a))| - |eqs(a)|  (failing >= -f + c):',f'{time.time()-t0:.0f}s')
for f,a,ne,ins in res[:25]:
    print(f'  a{a} fp={ne} |A|={len(ins)} f={f} A={list(ins)[:14]}')

# ---- greedy growth from the best seeds -----------------------------------
def grow(E0, nsteps=25):
    E=set(E0)
    # local pool: atoms sharing >=2 equations with E
    cnt=collections.Counter()
    for e in E:
        for b in eqat[e]: cnt[b]+=1
    pool={b for b,c in cnt.items() if c>=2}
    for _ in range(3):
        miss={b:EQS[b]-E for b in pool}
        addable=collections.Counter()
        for b,m in miss.items():
            if len(m)==1: addable[next(iter(m))]+=1
        if not addable: break
        e,g=addable.most_common(1)[0]
        if g-1<0: break
        E.add(e)
        for b in eqat[e]: pool.add(b)
    best=(len([b for b in pool if EQS[b]<=E])-len(E),frozenset(E))
    for _ in range(nsteps):
        cnt=collections.Counter()
        for e in E:
            for b in eqat[e]: cnt[b]+=1
        pool={b for b,c in cnt.items() if c>=2}
        addable=collections.Counter()
        for b in pool:
            m=EQS[b]-E
            if len(m)==1: addable[next(iter(m))]+=1
        if not addable: break
        e,g=addable.most_common(1)[0]
        E.add(e)
        A=[b for b in pool|set(eqat[e]) if EQS[b]<=E]
        f=len(A)-len(E)
        if f>best[0]: best=(f,frozenset(E))
        if g-1<0 and f<best[0]-2: break
    E=best[1]
    A=sorted({b for e in E for b in eqat[e] if EQS[b]<=E})
    return best[0],sorted(E),A

seen=set(); out=[]
for f,a,ne,ins in res[:400]:
    key=frozenset(EQS[a])
    if key in seen: continue
    seen.add(key)
    g,E,A=grow(EQS[a])
    out.append((g,a,len(E),len(A),E,A))
out.sort(key=lambda t:-t[0])
print('\nAFTER GROWTH (best f):')
for g,a,ne,na,E,A in out[:20]:
    print(f'  seed a{a}: f={g} |E|={ne} |S|={na}  ->  |E|-|S| = {ne-na}')
    print(f'      E={E}')
    print(f'      S={A}')
json.dump([[g,a,ne,na,E,A] for g,a,ne,na,E,A in out[:60]],open(os.path.join(HERE,'pa_scan2.json'),'w'))
print(f'{time.time()-t0:.0f}s done')
