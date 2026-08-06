"""CL step 1: exact linear forms mod p of the two cluster residues over ALL free inputs."""
import os, sys, json, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = 2**256 - 2**32 - 977

W = os.path.join(HERE, 'mod9118_0.json')
v = L.load(W)
vm = [x % P for x in v]
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
fail = L.failing_eqs(av)
print(f'frame mod9118_0: nonzero atoms {nz}, failing eqs {len(fail)}, score {L.NEQ-len(fail)}')
print(f'free inputs total: {len(ad.FREE)}')

CL1 = [21617, 29539]
forms = {}
for a in CL1:
    t0 = time.time()
    g = ad.grad(a, vm)
    forms[a] = g
    print(f'\na{a}: residue = {av[a]}  (mod p = {av[a]%P})')
    print(f'   gradient support {len(g)} free inputs  ({time.time()-t0:.1f}s)')
    print(f'   src: {L.atom_src[a][:200]}')
    print(f'   eqs: {sorted(L.atom2eq.get(a,{}))}')

# constant term:  residue(x) = const + sum g_u * x_u  (mod p)
for a in CL1:
    g = forms[a]
    s = sum(d*vm[u] for u,d in g.items()) % P
    const = (av[a] - s) % P
    print(f'a{a}: const term = {const}')
    forms[a]['__const__'] = const

json.dump({str(a): {str(u): str(d) for u,d in forms[a].items()} for a in CL1},
          open(os.path.join(HERE, 'cl_forms.json'),'w'))

# shared support?
s1 = set(forms[21617]) - {'__const__'}
s2 = set(forms[29539]) - {'__const__'}
print(f'\nsupport sizes: a21617 {len(s1)}, a29539 {len(s2)}, shared {len(s1&s2)}')
print(f'  shared: {sorted(s1&s2)[:40]}')
print(f'  only a21617: {sorted(s1-s2)[:40]}')
print(f'  only a29539: {len(s2-s1)}')

# also the second cluster
CL2 = [2423, 26731, 33929]
print('\n=== second cluster ===')
for a in CL2:
    g = ad.grad(a, vm)
    forms[a] = g
    print(f'a{a}: value {av[a]}  support {len(g)}  eqs {sorted(L.atom2eq.get(a,{}))}')
    print(f'   src: {L.atom_src[a][:200]}')
    sa = set(g)
    print(f'   overlap with a21617 support: {len(sa&s1)}, with a29539: {len(sa&s2)}')
json.dump({str(a): {str(u): str(d) for u,d in forms[a].items()} for a in forms},
          open(os.path.join(HERE, 'cl_forms.json'),'w'))
print('saved cl_forms.json')
