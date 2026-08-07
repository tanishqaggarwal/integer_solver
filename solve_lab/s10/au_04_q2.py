import os, sys, re, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256-2**32-977
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)

print('NA =', L.NA, ' atoms with >=1 equation:', len(L.atom2eq))
noeq = [a for a in range(L.NA) if a not in L.atom2eq]
print('atoms in ZERO equations:', len(noeq))
print('   of these, gate outputs (definers):', sum(1 for a in noeq if a in L.atom_out))
print('   of these, checks:', sum(1 for a in noeq if a not in L.atom_out))
print('   nonzero-valued among them:', sum(1 for a in noeq if av[a]))

# --- eq 8680 (the only equation containing a37887) ---
m,sq,co = L.eq_atoms[8680]
print(f'\n=== eq 8680: mult={m} sq={sq} n_atoms={len(co)}')
for a,c in sorted(co.items()):
    ne = len(L.atom2eq.get(a,{}))
    print(f'   a{a:<6} coeff={c:<5} n_eqs={ne:<4} nonzero={av[a]!=0} gate_out={L.atom_out.get(a)}')

# --- parse Q's constituent atoms by source matching ---
src2a = collections.defaultdict(list)
for a in range(L.NA):
    src2a[L.atom_src[a].strip()].append(a)

s = L.atom_src[37887]
# split the squared source into the single factor
half = s[1:len(s)//2]   # crude; instead find the ') * (' split
i = s.find(') * (')
f1 = s[1:i]
print('\n=== Q factor length', len(f1))
# top-level split of f1 on ' + ' with paren depth 0
parts=[]; depth=0; cur=''
k=0
while k < len(f1):
    ch=f1[k]
    if ch=='(': depth+=1
    if ch==')': depth-=1
    if depth==0 and f1[k:k+3]==' + ':
        parts.append(cur); cur=''; k+=3; continue
    cur+=ch; k+=1
parts.append(cur)
print(f'{len(parts)} top-level terms in Q')
tot=[]
for pt in parts:
    pt=pt.strip()
    mm = re.match(r'^(-?\d+) \* \((.*)\)$', pt)
    if mm: coef=int(mm.group(1)); body=mm.group(2)
    else: coef=1; body=pt
    cands = src2a.get(body.strip(), [])
    ne = [len(L.atom2eq.get(a,{})) for a in cands]
    tot.append((coef, body[:60], cands, ne))
    print(f'  coeff {coef:>4}  atoms={cands}  n_eqs={ne}  nonzero={[av[a]!=0 for a in cands]}  src={body[:70]}')
