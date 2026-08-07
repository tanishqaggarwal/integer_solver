"""Census of alternative orientations.

`fwd2.target_of` takes the LEADING bare variable of `x_A - rest` as the atom's target.  But an atom
is solvable for ANY variable that occurs as a top-level additive term with coefficient +-1 -- the
resulting definition is still integral, so it is an equally legal orientation.  This enumerates
those alternatives per atom, which is the raw material for re-orienting the frame.

Top-level additive decomposition: walk the Add/Sub spine of the expression and keep the terms that
are a bare Name (sign +-1).  Anything under a Mult, or with a numeric coefficient, is not a legal
unit target and is ignored.
"""
import ast, json, pickle, os, sys
from collections import Counter, defaultdict
import model

HERE = os.path.dirname(os.path.abspath(__file__))
d = model.get()
atom_src = d['atom_src']
atom_vars = d['atom_vars']
eq_terms = d['eq_terms']
NA = len(atom_src)
NV = 38748


def spine(node, sign, out):
    """collect (sign, node) over the top-level +/- spine"""
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        spine(node.left, sign, out)
        spine(node.right, sign if isinstance(node.op, ast.Add) else -sign, out)
    elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        spine(node.operand, -sign, out)
    else:
        out.append((sign, node))


def unit_targets(src):
    """variables occurring as a bare top-level term with coefficient +-1 -> {var: sign}"""
    t = ast.parse(src, mode='eval').body
    out = []
    spine(t, 1, out)
    res = {}
    for s, n in out:
        if isinstance(n, ast.Name) and n.id.startswith('x_'):
            v = int(n.id[2:])
            if v in res:            # appears twice at top level: not a clean unit target
                res[v] = None
            else:
                res[v] = s
    return {k: v for k, v in res.items() if v is not None}


if __name__ == '__main__':
    alts = {}
    cnt = Counter()
    for a in range(NA):
        u = unit_targets(atom_src[a])
        alts[a] = u
        cnt[len(u)] += 1
    print('atoms: %d' % NA)
    print('number of legal unit targets per atom: %s' % sorted(cnt.items()))
    multi = [a for a in range(NA) if len(alts[a]) > 1]
    print('atoms admitting MORE THAN ONE orientation: %d (%.1f%%)' % (len(multi), 100 * len(multi) / NA))

    F = pickle.load(open(os.path.join(HERE, 'fwd2.pkl'), 'rb'))
    tgt = F['tgt']
    definer = F['definer']
    checks = set(F['checks'])
    free0 = set(F['free0'])
    print('\ncurrent frame: %d definitions, %d check atoms, %d pure free inputs'
          % (sum(1 for v in range(NV) if definer[v] >= 0), len(checks), len(free0)))

    used = [a for a in range(NA) if a not in checks]
    print('definition atoms with an alternative target: %d of %d'
          % (sum(1 for a in used if len(alts[a]) > 1), len(used)))

    REGION_ATOMS = [22229, 22230, 22231, 35758, 35759, 35760, 35761, 35762, 37887]
    print('\nregion atoms and their legal orientations:')
    for a in REGION_ATOMS:
        print('  atom %-6d cur_tgt=%-7s check=%-5s  unit targets: %s'
              % (a, tgt[a], a in checks, alts[a]))

    POOL = json.load(open(os.path.join(HERE, 'pool.json')))
    KEY = [642, 28730, 29854, 31864]
    print('\nthe 4 live pool variables, their definers and that definer\'s alternatives:')
    for v in KEY:
        a = definer[v]
        print('  x_%-6d definer atom %-6d  alternatives: %s' % (v, a, alts[a] if a >= 0 else None))
    npool_alt = sum(1 for v in POOL if definer[v] >= 0 and len(alts[definer[v]]) > 1)
    print('\npool variables whose definer admits an alternative target: %d of %d'
          % (npool_alt, len(POOL)))
    json.dump({str(k): {str(kk): vv for kk, vv in v.items()} for k, v in alts.items()},
              open(os.path.join(HERE, 'alts.json'), 'w'))
    print('wrote alts.json')
