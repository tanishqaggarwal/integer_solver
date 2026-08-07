"""Efficient global deficiency scan: for every atom a, count atoms fully inside eqs(a),
plus the cheapest knob-controlled atoms."""
import os, sys, collections, json, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
EQS=[frozenset(L.atom2eq[a]) for a in range(L.NA)]
eqat=[sorted(co) for m,sq,co in L.eq_atoms]
free=set(u for u in range(L.NVARS) if u not in L.definer)
knobatom={}
for u in free:
    if len(L.var_atoms[u])==1: knobatom.setdefault(L.var_atoms[u][0],[]).append(u)
print('knob atoms',len(knobatom))
cheap=sorted(knobatom,key=lambda a:len(EQS[a]))[:20]
print('cheapest knob-controlled atoms:')
for a in cheap:
    print(f'  a{a} fp={len(EQS[a])} knobs={knobatom[a]} eqs={sorted(EQS[a])}  {L.atom_src[a][:80]}')
t0=time.time()
best=[]
for a in range(L.NA):
    Ea=EQS[a]
    cnt=collections.Counter()
    for e in Ea:
        for b in eqat[e]: cnt[b]+=1
    ins=[b for b,c in cnt.items() if c==len(EQS[b])]
    f=len(ins)-len(Ea)
    best.append((f,a,len(Ea),len(ins),sorted(ins)))
    if a%8000==0: print(a,f'{time.time()-t0:.0f}s',flush=True)
best.sort(key=lambda t:-t[0])
print('\nTOP single-seed deficiency (f = |A(eqs(a))| - |eqs(a)|):')
for f,a,ne,ni,ins in best[:30]:
    print(f'  a{a} fp={ne} atoms_inside={ni} f={f} inside={ins[:20]}')
json.dump([[f,a,ne,ni,ins] for f,a,ne,ni,ins in best[:600]],open(os.path.join(HERE,'pa_scan.json'),'w'))
print(f'{time.time()-t0:.0f}s done')
