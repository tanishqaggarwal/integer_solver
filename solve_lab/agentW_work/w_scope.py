"""W stage 13: scope boundary.  (i) do the 766 off-pins also hold in the deliverable?
(ii) what are the 7 failing atoms?  (iii) how exposed are the 1149 congruence atoms to
equation-level cancellation (the one freedom the atom-level theorem does NOT cover)?"""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB, model
from collections import Counter
d = model.get(); A = d['atom_src']; AV = d['atom_vars']; EQ = d['eq_terms']
byvar = {}
for i, vs in enumerate(AV):
    for v in vs: byvar.setdefault(v, []).append(i)
blocks = json.load(open('w_blocks4.json'))
fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st = frameB.State(fr, fv); v = fr.forward(fv)
print('score', st.score())
def short(u): return [a for a in byvar.get(u, []) if len(A[a]) < 200]
PATS = [re.compile(r'x_(\d+) \* x_(\d+) - x_(\d+)$'), re.compile(r'x_(\d+) \* x_(\d+) \+ x_(\d+)$'),
        re.compile(r'(-?\d+) \* \(x_(\d+) \* x_(\d+)\) - x_(\d+)$'), re.compile(r'x_(\d+) \* x_(\d+) - (-?\d+) \* x_(\d+)$')]
nz = 0; tot = 0
for b in blocks:
    L = b['L']
    NL = [int(re.match(r'x_(\d+)', A[a]).group(1)) for a in short(L)
          if re.fullmatch(r'x_\d+ - \(1 - x_%d\)' % L, A[a])][0]
    for iv in (b['i5'], b['i6']):
        for a in short(iv):
            s = A[a]
            if not any(p.fullmatch(s) for p in PATS): continue
            vs = set(int(x) for x in re.findall(r'x_(\d+)', s))
            if not ({NL, iv} <= vs): continue
            tot += 1
            if eval(frameB.ACODE[a], {'v': v, '__builtins__': {}}) != 0: nz += 1
print('off-pin atoms evaluated in the deliverable: %d, nonzero: %d' % (tot, nz))
# the 7 failing atoms
print()
NZ = sorted(st.nz())
for a in NZ:
    print('  a%-6d %s' % (a, A[a][:170]))
# incidence
eq_of = Counter()
for i, (m, sq, tl) in enumerate(EQ):
    for c, a in tl: eq_of[a] += 1
cong = [cg['ring']['atom'] for b in blocks for cg in b['congs']]
print()
print('congruence atoms: %d ; #equations each appears in:' % len(cong), Counter(eq_of[a] for a in cong).most_common(8))
print('is a congruence atom ever ALONE in an equation?',
      sum(1 for i,(m,sq,tl) in enumerate(EQ) if len(tl)==1 and tl[0][1] in set(cong)))
print('mean eqs per congruence atom: %.2f' % (sum(eq_of[a] for a in cong)/len(cong)))
