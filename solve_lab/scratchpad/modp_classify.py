import sys, os, json
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/scratchpad')
import atomlib as A
p = A.p

# build wire set (signed copy-class of x_26064) - reuse lin_extract logic
atoms = A.ATOMS
par = list(range(A.NVARS));
def find(x):
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
def uni(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: par[rb] = ra
for poly in atoms:
    if all(len(vs) <= 1 for vs, c in poly):
        lin = [(vs[0], c) for vs, c in poly if len(vs) == 1]
        const = sum(c for vs, c in poly if len(vs) == 0)
        if len(lin) == 2 and const == 0 and abs(lin[0][1]) == 1 and abs(lin[1][1]) == 1:
            uni(lin[0][0], lin[1][0])
r = find(26064)
wire = set(v for v in range(A.NVARS) if find(v) == r)
print(f"wire members: {len(wire)}")

# classify each atom mod p: drop wire linear terms (=0) and wire-containing product terms (=0)
n_lin = n_quad = n_const_only = n_higher = 0
quad_atoms = []      # genuinely quadratic mod p
quad_monomials = set()
nonwire_vars_in_quad = set()
for ai, poly in enumerate(atoms):
    max_deg = 0
    has_nonwire_prod = False
    prod_pairs = []
    for vs, c in poly:
        # reduce: a term with any wire var -> the term is a product; if it has a wire factor and it's a product (len>=2) it's 0 mod p (wire=0). If len==1 wire -> 0.
        if len(vs) == 0:
            continue
        if any(v in wire for v in vs):
            # term vanishes mod p (wire ~ 0)  (whether linear wire or product-with-wire)
            continue
        # non-wire term
        if len(vs) == 1:
            max_deg = max(max_deg, 1)
        elif len(vs) == 2:
            max_deg = max(max_deg, 2)
            has_nonwire_prod = True
            prod_pairs.append(tuple(sorted(vs)))
        else:
            max_deg = max(max_deg, len(vs))
    if max_deg >= 3:
        n_higher += 1
    elif has_nonwire_prod:
        n_quad += 1
        quad_atoms.append(ai)
        for pp in prod_pairs:
            quad_monomials.add(pp); nonwire_vars_in_quad.update(pp)
    elif max_deg == 1:
        n_lin += 1
    else:
        n_const_only += 1

print(f"\nAtom classification mod p (wire terms dropped):")
print(f"  linear:            {n_lin}")
print(f"  genuinely quadratic:{n_quad}")
print(f"  higher (>=3):      {n_higher}")
print(f"  const-only/zero:   {n_const_only}")
print(f"\nGenuine quadratic monomials (distinct products): {len(quad_monomials)}")
print(f"Distinct non-wire vars appearing in a genuine product: {len(nonwire_vars_in_quad)}")
# show a few quad atoms
print(f"\nSample genuine-quadratic atoms:")
for ai in quad_atoms[:15]:
    print(f"  atom {ai}: {A.ATOM_REPR[ai][:90]}")
json.dump({'wire':sorted(wire),'quad_atoms':quad_atoms,
           'quad_monomials':[list(m) for m in quad_monomials],
           'nonwire_vars_in_quad':sorted(nonwire_vars_in_quad)},
          open('/home/user/integer_solver/solve_lab/scratchpad/modp_class.json','w'))
