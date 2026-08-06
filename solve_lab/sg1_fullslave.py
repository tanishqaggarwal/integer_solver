"""Iterative full-slave forward: forward() then re-slave every checked free input whose check-eq
residual != 0, to fixpoint. Keeps all equality-checks satisfied under perturbation of unchecked frees."""
import sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p
from collections import defaultdict, deque

NV = H.NVARS; freeinp = H.freeinp; anc = H.anc
eqvars = H.eqvars; eqcode = H.eqcode; order = H.order; val = H.val

sm = json.load(open('sg1_slavemap.json'))
slave_of = {int(v): (i, c) for v,(i,c) in sm.items()}
checked = set(slave_of)
ns = {'v': val, '__builtins__': {}}

# order checked frees by depth = number of gate-ancestors of their check-eq (heuristic: shallow first)
def eq_depth(i):
    return len(eqvars[i])
slave_order = sorted(slave_of, key=lambda v: eq_depth(slave_of[v][0]))

def forward_iter_slave(max_iter=60):
    baddiv = 0
    for it in range(max_iter):
        H.forward()
        changed = 0; baddiv = 0
        for v in slave_order:
            i, c = slave_of[v]
            res = eval(eqcode[i], ns)
            if res == 0: continue
            if res % c == 0:
                val[v] -= res // c; changed += 1
            else:
                baddiv += 1
        if changed == 0:
            H.forward()
            return it, baddiv
    H.forward()
    return max_iter, baddiv

def fails():
    return [i for i,co in enumerate(eqcode) if eval(co, ns)!=0]

if __name__ == '__main__':
    v013 = H.loadd('best/new_instance_partial_39013.json')
    for v in freeinp: val[v] = v013.get(v, 0)
    it, bd = forward_iter_slave()
    F = fails()
    print(f"39013 iter-slave: {len(F)} fail (converged in {it} iters, baddiv={bd})")
    print(f"  fails: {sorted(F)}")

    # residue fix
    r29 = val[29322] % p; r35 = val[3558] % p
    val[14853] -= r29; val[16742] += r35
    it2, bd2 = forward_iter_slave()
    F2 = fails()
    print(f"\nafter residue fix + iter-slave: {len(F2)} fail (iters={it2}, baddiv={bd2})")
    print(f"  fails: {sorted(F2)}")
    print(f"  x_29322%p={val[29322]%p}, x_3558%p={val[3558]%p}")
    print(f"  G1={7376877*val[642]+val[2099]-val[7068]}, G2={val[4432]-val[19964]-val[28730]}")
