"""W stage 14: exactly WHICH atoms does the deliverable break, in the block vocabulary?"""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB, model
from collections import Counter
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
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
def short(u): return [a for a in byvar.get(u, []) if len(A[a]) < 200]
PATS = [re.compile(r'x_(\d+) \* x_(\d+) - x_(\d+)$'), re.compile(r'x_(\d+) \* x_(\d+) \+ x_(\d+)$'),
        re.compile(r'(-?\d+) \* \(x_(\d+) \* x_(\d+)\) - x_(\d+)$'), re.compile(r'x_(\d+) \* x_(\d+) - (-?\d+) \* x_(\d+)$')]
OFF = {}
for b in blocks:
    L = b['L']
    NL = [int(re.match(r'x_(\d+)', A[a]).group(1)) for a in short(L)
          if re.fullmatch(r'x_\d+ - \(1 - x_%d\)' % L, A[a])][0]
    for slot, iv in (('i5', b['i5']), ('i6', b['i6'])):
        for a in short(iv):
            s = A[a]
            if not any(p.fullmatch(s) for p in PATS): continue
            if not ({NL, iv} <= set(int(x) for x in re.findall(r'x_(\d+)', s))): continue
            OFF[a] = (b['E'], slot, iv, L, NL)
CONG = {cg['ring']['atom']: (b['E'], 'cong') for b in blocks for cg in b['congs']}
# leaf pins:  sel*(w - K) - m*z
LEAF = {}
for i, s in enumerate(A):
    m = re.fullmatch(r'x_(\d+) \* \(x_(\d+) - (\d+)\) - (?:(\d+) \* )?x_(\d+)', s)
    if m and len(m.group(3)) > 60: LEAF[i] = (int(m.group(1)), int(m.group(2)))
print('atom inventory: %d off-pins, %d congruences, %d leaf pins' % (len(OFF), len(CONG), len(LEAF)))
print()
NZ = sorted(st.nz())
for a in NZ:
    tag = 'OFF-PIN %s' % (OFF[a],) if a in OFF else ('CONGRUENCE %s' % (CONG[a],) if a in CONG
          else ('LEAF-PIN %s' % (LEAF[a],) if a in LEAF else 'other'))
    val = eval(frameB.ACODE[a], {'v': v, '__builtins__': {}})
    print('  a%-6d  %-28s  val%%P==0: %s   %s' % (a, tag, val % PVAL == 0, A[a][:80]))
print()
nzoff = [a for a in OFF if eval(frameB.ACODE[a], {'v': v, '__builtins__': {}}) != 0]
print('nonzero off-pins:', [(a, OFF[a]) for a in nzoff])
for a in nzoff:
    E, slot, iv, L, NL = OFF[a]
    print('   block E=%d slot %s: gate L=x_%d -> %s ,  NOT=x_%d -> %s ,  i=x_%d -> %s (mod P: %s)'
          % (E, slot, L, v[L], NL, v[NL], iv, str(v[iv])[:40], v[iv] % PVAL == 0))
nzc = [a for a in CONG if eval(frameB.ACODE[a], {'v': v, '__builtins__': {}}) != 0]
print('nonzero congruence atoms:', nzc)
nzl = [a for a in LEAF if eval(frameB.ACODE[a], {'v': v, '__builtins__': {}}) != 0]
print('nonzero leaf pins:', nzl)
