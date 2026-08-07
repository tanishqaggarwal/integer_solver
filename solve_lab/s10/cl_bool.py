"""CL: the cone of x_27522 is a boolean mux tree.  Enumerate the boolean strata."""
import os, sys, json, collections, itertools, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import cl_engine as E
P = E.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE)

v0 = L.load(os.path.join(HERE,'mod9118_0.json'))
av0, nz0, S0, bad0 = E.stats(v0)
print(f'base score {S0} nz {nz0}')

def cone(u):
    seen=set(); st=[u]
    while st:
        w=st.pop()
        if w in seen: continue
        seen.add(w)
        a = definer.get(w)
        if a is None: continue
        for z in L.avars[a]:
            if z!=w and z not in seen: st.append(z)
    return seen

C1 = cone(27522)
F1 = sorted(u for u in C1 if u in FREE)
print(f'\ncone(x_27522): {len(C1)} vars, {len(F1)} free inputs')
for u in F1:
    print(f'  x_{u:<6} val={str(v0[u])[:30]:<32} bits={v0[u].bit_length():<5} consumers={len(L.var_atoms[u])}')

BOOLS = [u for u in F1 if v0[u] in (0,1)]
print(f'\nboolean-valued free inputs in cone: {len(BOOLS)} -> {BOOLS}')

# the mux selectors
SEL = [28940, 23047]
PAY = [19799, 36462, 8239]
print(f'\ncurrent selectors: ' + ', '.join(f'x_{s}={v0[s]}' for s in SEL))
print(f'payloads: ' + ', '.join(f'x_{p}={"FREE " if p in FREE else ""}{str(v0[p])[:20]}' for p in PAY))

# forward-evaluate JUST the cone for a given boolean setting
order1 = [t for t in ad.ORDER if t in C1 and t in definer]
def eval_cone(setting):
    v = list(v0)
    for u, b in setting.items(): v[u] = b
    for t in order1:
        nv = T.solve_lin(definer[t], t, v)
        if nv is not None: v[t] = nv
    return v

print('\n=== enumerating boolean strata of cone(x_27522) ===')
res = collections.defaultdict(list)
n = len(BOOLS)
if n <= 14:
    for bits in itertools.product((0,1), repeat=n):
        st = dict(zip(BOOLS, bits))
        v = eval_cone(st)
        key = (v[28940], v[23047])
        res[key].append(bits)
    for k, lst in sorted(res.items()):
        print(f'  selectors {k}: {len(lst)} settings; example {dict(zip(BOOLS,lst[0]))}')
        # what would x_27522 be
        v = eval_cone(dict(zip(BOOLS, lst[0])))
        print(f'      x_27522 = {str(v[27522])[:40]}  ({v[27522].bit_length()} bits)  mod p = {v[27522]%P}')
        print(f'      target x_14623 mod p = {v0[14623]%P}')
json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE,'cl_bool.json'),'w'))
