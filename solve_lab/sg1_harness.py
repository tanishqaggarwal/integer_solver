"""sg1 solver harness: slaved forward + residue fix experiments.
Keeps x_7068:=x_2099 (G1) and x_4432:=x_19964 (G2) slaved. Never moves them independently."""
import sys, json
sys.path.insert(0, '/home/user/integer_solver/solve_lab')
import heal_harness as H
p = H.p
from collections import defaultdict

# slaved pairs: free_input := gate_output  (keeps the equality-check atom satisfied)
SLAVE = {7068: 2099, 4432: 19964}

def forward_slaved(iters=3):
    """forward(), then set slaved free inputs to their RHS gate outputs, re-forward. Repeat."""
    for _ in range(iters):
        H.forward()
        changed = False
        for f, g in SLAVE.items():
            if H.val[f] != H.val[g]:
                H.val[f] = H.val[g]; changed = True
        if not changed:
            break
    H.forward()

def load_state(path):
    d = H.loadd(path)
    for v in H.freeinp: H.val[v] = d.get(v, 0)
    forward_slaved()

def fails():
    return H.fails()

def save(path, note=""):
    d = {f"x_{i}": H.val[i] for i in range(H.NVARS)}
    json.dump(d, open(path, 'w'))
    if note: print(f"saved {path}: {note}")

# eq free-support precompute
_eq_free_support = None
def eq_free_support():
    global _eq_free_support
    if _eq_free_support is None:
        _eq_free_support = []
        for i in range(len(H.eqvars)):
            s = set()
            for var in H.eqvars[i]:
                if var in H.freeinp: s.add(var)
                else: s |= H.anc.get(var, set())
            _eq_free_support.append(s)
    return _eq_free_support

if __name__ == '__main__':
    load_state('best/new_instance_partial_39013.json')
    F0 = set(fails())
    print(f"39013 slaved baseline: {len(F0)} fail: {sorted(F0)}")
    # G1/G2 check
    G1 = 7376877*H.val[642] + H.val[2099] - H.val[7068]
    G2 = H.val[4432] - H.val[19964] - H.val[28730]
    print(f"G1={G1}, G2={G2}")
    print(f"x_29322%p = {H.val[29322]%p}")
    print(f"x_3558%p  = {H.val[3558]%p}")
