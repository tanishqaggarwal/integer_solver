import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977

SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
fail = L.failing_eqs(av)
print('score', L.NEQ - len(fail), 'failing', fail)

# --- (b) recompute E from scratch ---
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
print('E from atom2eq:', E, len(E))
# independent re-derivation directly from eq_terms, summing duplicates
E2 = []
for i,(m,sq,tl) in enumerate(L.eq_terms):
    co = collections.defaultdict(int)
    for c,a in tl: co[a]+=c
    for a in SEVEN:
        if co.get(a,0) != 0:
            E2.append(i); break
print('E2 (nonzero coeff after dedup):', E2, len(E2))
# equations where a SEVEN atom appears syntactically but cancels to zero
Ecancel=[]
for i,(m,sq,tl) in enumerate(L.eq_terms):
    co = collections.defaultdict(int)
    present=set()
    for c,a in tl:
        co[a]+=c
        if a in SEVEN: present.add(a)
    if present and not any(co.get(a,0) for a in present):
        Ecancel.append(i)
print('equations where a SEVEN atom cancels to zero coeff:', Ecancel)

# --- multiplier / is_square check on those 12 ---
print('\n--- the 12 equations in full ---')
allatoms = set()
for e in E:
    m,sq,co = L.eq_atoms[e]
    S = sum(c*av[a] for a,c in co.items())
    val = m*(S*S if sq else S)
    others = {a:c for a,c in co.items() if a not in SEVEN}
    nzother = {a:(c,av[a]) for a,c in others.items() if av[a]}
    allatoms |= set(co)
    print(f' eq {e:>6} mult={m} sq={sq} natoms={len(co)} SEVENrow={[co.get(a,0) for a in SEVEN]}'
          f' n_other={len(others)} nzother={ {k:(c,str(x)[:12]) for k,(c,x) in nzother.items()} } S!=0:{S!=0} val!=0:{val!=0}')
print('\nmults of all 39033 eqs: distinct =', sorted(set(m for m,sq,tl in L.eq_terms))[:20])
print('any mult == 0 anywhere?', any(m==0 for m,sq,tl in L.eq_terms))
print('n is_square eqs:', sum(1 for m,sq,tl in L.eq_terms if sq))
print('is_square among the 12:', [(e, L.eq_atoms[e][1]) for e in E])

print('\nall atoms appearing in the 12 equations:', len(allatoms))
# for each such atom: how many equations OUTSIDE the 12 does it touch, and its value
info=[]
for a in sorted(allatoms):
    eqs = set(L.atom2eq.get(a,{}))
    out = len(eqs - set(E))
    inn = len(eqs & set(E))
    info.append((out, inn, a, av[a]!=0, a in L.atom_out))
info.sort()
for out,inn,a,nz,isout in info:
    print(f'  a{a:<6} in12={inn:<3} outside={out:<4} nonzero={nz} is_gate_out={isout}')
