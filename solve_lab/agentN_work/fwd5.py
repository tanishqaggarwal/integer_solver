"""Rebuild the frame wholesale under a different target rule.

`fwd2` orients every atom by its LEADING bare variable.  An atom `x_v - rest` is legally solvable
for ANY variable occurring as a top-level additive term with coefficient +-1, so there are many
global orientations.  This builds them.

Definition semantics without re-parsing: an atom's value is  s*x_v + F(other vars)  with s = +-1
for a legal unit target v, so setting the atom to 0 means  x_v <- x_v - s*atom_value.  Evaluating
the atom with every other variable already known therefore realises the definition exactly, and a
topological order is produced by the same greedy propagation `fwd2` uses.

A frame is judged by what it forces to zero: definition atoms are identically zero for every choice
of free inputs, so only CHECK atoms carry value into the equations.  The question is whether some
orientation leaves fewer check atoms carrying nonzero value in the failing equations -- in
particular whether the 7 nonzero atoms {22229,22230,35758,35759,35760,35761,35762} survive.
"""
import ast, json, os, pickle, random, re, sys, time
from collections import defaultdict, deque
import model
from orient import unit_targets

HERE = os.path.dirname(os.path.abspath(__file__))
VAR_RE = re.compile(r'x_(\d+)')
d = model.get()
atom_src = d['atom_src']
atom_vars = d['atom_vars']
eq_terms = d['eq_terms']
NA = len(atom_src)
NV = 38748

print('computing legal unit targets for %d atoms ...' % NA, flush=True)
_t = time.time()
U = [unit_targets(atom_src[a]) for a in range(NA)]
print('  %.0fs' % (time.time() - _t), flush=True)

ATOMCODE = {}


def code(a):
    c = ATOMCODE.get(a)
    if c is None:
        c = compile(VAR_RE.sub(r'v[\1]', atom_src[a]), '<a>', 'eval')
        ATOMCODE[a] = c
    return c


W = json.load(open(os.path.join(HERE, '..', 'best', 'new_instance_partial_39026.json')))
wv = [0] * NV
for k, val in W.items():
    wv[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)

var_atoms = defaultdict(list)
for a in range(NA):
    for v in atom_vars[a]:
        var_atoms[v].append(a)

# a variable can only ever be defined if it is a legal unit target of some atom
targetable = set()
for a in range(NA):
    targetable |= set(U[a])
PURE_FREE = [v for v in range(NV) if v not in targetable]
print('variables that are a legal unit target of some atom: %d ; never targetable: %d'
      % (len(targetable), len(PURE_FREE)), flush=True)


def choose_targets(rule, seed=0, prefer=None):
    """One chosen target per atom, exactly as fwd2 fixes tgt[] before propagating."""
    rnd = random.Random(seed)
    tgt = [None] * NA
    sgn = [0] * NA
    for a in range(NA):
        ks = list(U[a].items())
        if not ks:
            continue
        if rule == 'first':
            pick = ks[0]
        elif rule == 'last':
            pick = ks[-1]
        elif rule == 'random':
            pick = rnd.choice(ks)
        elif rule == 'lowvar':
            pick = min(ks, key=lambda z: z[0])
        elif rule == 'highvar':
            pick = max(ks, key=lambda z: z[0])
        elif rule == 'prefer':
            cand = [z for z in ks if prefer and z[0] in prefer]
            pick = cand[0] if cand else ks[0]
        else:
            raise ValueError(rule)
        tgt[a], sgn[a] = pick
    return tgt, sgn


def build(rule, seed=0, prefer=None):
    """fwd2's propagation, but over a chosen target map.
    Returns (defs, checks, free0), defs a topological list of (var, atom, sign)."""
    tgt, sgn = choose_targets(rule, seed, prefer)
    istarget = set(t for t in tgt if t is not None)
    known = bytearray(NV)
    for v in range(NV):
        if v not in istarget:
            known[v] = 1
    unk = [0] * NA
    for a in range(NA):
        t = tgt[a]
        unk[a] = sum(1 for v in atom_vars[a] if v != t and not known[v])
    definer = [-1] * NV
    order = []
    Q = deque(a for a in range(NA) if tgt[a] is not None and unk[a] == 0)
    while Q:
        a = Q.popleft()
        t = tgt[a]
        if t is None or known[t] or unk[a] != 0:
            continue
        known[t] = 1
        definer[t] = a
        order.append((t, a, sgn[a]))
        for b in var_atoms[t]:
            if tgt[b] != t:
                unk[b] -= 1
                if unk[b] == 0 and tgt[b] is not None and not known[tgt[b]]:
                    Q.append(b)
    used = set(definer[v] for v in range(NV) if definer[v] >= 0)
    checks = [a for a in range(NA) if a not in used]
    free0 = [v for v in range(NV) if definer[v] < 0]
    return order, checks, free0


def evaluate(order, checks, free0, freevals=None):
    """forward-propagate from the witness values on the free inputs"""
    v = list(wv) if freevals is None else list(freevals)
    ns = {'v': v, '__builtins__': {}}
    for var, a, s in order:
        v[var] = v[var] - s * eval(code(a), ns)
    av = {}
    for a in checks:
        av[a] = eval(code(a), ns)
    fails = []
    for i, (m, sq, tl) in enumerate(eq_terms):
        t = 0
        for c, a in tl:
            x = av.get(a)
            if x:
                t += c * x
        if (m * (t * t if sq else t)) != 0:
            fails.append(i)
    return v, av, fails


REGION7 = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
REGION9 = REGION7 + [22231, 37887]

if __name__ == '__main__':
    rules = [('first', 0), ('last', 0), ('lowvar', 0), ('highvar', 0),
             ('random', 1), ('random', 2), ('random', 3), ('random', 4), ('random', 5)]
    # a rule that deliberately tries to make the region atoms definitions
    prefer = set()
    for a in REGION9:
        prefer |= set(U[a])
    rules.append(('prefer', 0))
    out = []
    for rule, seed in rules:
        t0 = time.time()
        order, checks, free0 = build(rule, seed, prefer)
        tb = time.time() - t0
        t0 = time.time()
        v, av, fails = evaluate(order, checks, free0)
        te = time.time() - t0
        regdef = [a for a in REGION9 if a not in set(checks)]
        regchk = [a for a in REGION9 if a in set(checks)]
        nzchecks = [a for a in checks if av.get(a, 0) != 0]
        nz_in_region = [a for a in REGION9 if av.get(a, 0) != 0]
        tag = '%s%s' % (rule, ('/%d' % seed) if rule == 'random' else '')
        print('\n=== rule %-10s defs=%-6d checks=%-6d free=%-5d  score=%d  (build %.0fs eval %.0fs)'
              % (tag, len(order), len(checks), len(free0), len(eq_terms) - len(fails), tb, te),
              flush=True)
        print('   failing equations: %d %s' % (len(fails), fails[:12]), flush=True)
        print('   nonzero CHECK atoms anywhere: %d' % len(nzchecks), flush=True)
        print('   of the 9 region atoms: %d are definitions %s' % (len(regdef), regdef), flush=True)
        print('   region atoms still CHECKS: %s' % regchk, flush=True)
        print('   region atoms NONZERO here: %s' % nz_in_region, flush=True)
        out.append(dict(rule=tag, defs=len(order), checks=len(checks), free=len(free0),
                        score=len(eq_terms) - len(fails), nfail=len(fails), fails=fails[:40],
                        region_defs=regdef, region_checks=regchk, region_nonzero=nz_in_region,
                        nz_checks=len(nzchecks)))
        json.dump(out, open(os.path.join(HERE, 'runs', 'fwd5.json'), 'w'), indent=1)
    print('\nbest forward score over the orientations tried: %d' % max(r['score'] for r in out))
