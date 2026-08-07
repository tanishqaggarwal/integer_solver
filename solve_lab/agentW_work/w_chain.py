"""W stage 16: what actually connects the 7 broken atoms to the degenerate block E=33469?
Trace 33469's four inputs back through the definition DAG to the corrupted variables."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB, model
from collections import deque
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st = frameB.State(fr, fv); v = fr.forward(fv)
blocks = {b['E']: b for b in json.load(open('w_blocks4.json'))}
b = blocks[33469]
print('degenerate block E=33469  gate L=x_%d -> %s' % (b['L'], v[b['L']]))
for k in ('i1','i2','i3','i4','i5','i6'):
    print('   %s = x_%-6d  val%%P = %s' % (k, b[k], str(v[b[k]] % PVAL)[:44]))
print('   A mod P =', (v[b['i1']]-v[b['i2']]) % PVAL, '  B mod P =', (v[b['i4']]-v[b['i3']]) % PVAL)
# ancestors of the four live inputs, restricted to the definition DAG
definer = fr.definer
DET = {642, 28730, 29854, 31864}
def anc(v0_, cap=200000):
    seen = set(); Q = deque([v0_])
    while Q:
        x = Q.popleft()
        if x in seen or len(seen) > cap: continue
        seen.add(x)
        a = definer[x]
        if a < 0: continue
        for u in AV[a]:
            if u != x and u not in seen: Q.append(u)
    return seen
for k in ('i1','i2','i3','i4'):
    s = anc(b[k])
    print('   %s = x_%-6d : ancestors %6d ; contains detached %s' % (k, b[k], len(s), sorted(s & DET)))
# which block produces 33469's X and Y slots?  find blocks whose mux output feeds them
print()
print('the deliverable perturbs 4 free vars; what do they reach?')
for u in sorted(DET):
    print('   x_%-6d free=%s  value=%s' % (u, u in fr.free, str(v[u])[:50]))
    print('        descendants: %d   check atoms it can move: %d' % (len(fr.desc.get(u, [])), len(fr.chk.get(u, []))))
    print('        check atoms:', sorted(fr.chk.get(u, []))[:12])
