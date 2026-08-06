"""Session-9 harness: raw-equation probing around the best partial."""
import sys, re, json, time, os, pickle
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EQ_PATH = os.path.join(BASE, '..', 'EQUATIONS.txt')
NVARS = 38748
VAR_RE = re.compile(r'x_(\d+)')
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eqcache.pkl')

def load_raw():
    lines = []
    with open(EQ_PATH) as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line.rsplit('=', 1)[0])
    return lines

def load_equations():
    lines = load_raw()
    codes = [compile(VAR_RE.sub(r'v[\1]', s), '<eq>', 'eval') for s in lines]
    varsets = [frozenset(int(m) for m in VAR_RE.findall(s)) for s in lines]
    return codes, varsets

def load_assignment(path):
    with open(path) as f:
        d = json.load(f)
    v = [0]*NVARS
    for k, val in d.items():
        idx = int(k[2:]) if k.startswith('x_') else int(k)
        v[idx] = int(val)
    return v

def save_assignment(v, path):
    json.dump({f"x_{i}": v[i] for i in range(NVARS) if v[i] != 0}, open(path,'w'))

def evaluate(codes, v, idxs=None):
    ns = {'v': v, '__builtins__': {}}
    if idxs is None:
        return [i for i,c in enumerate(codes) if eval(c, ns) != 0]
    return [i for i in idxs if eval(codes[i], ns) != 0]

def resid(codes, v, i):
    return eval(codes[i], {'v': v, '__builtins__': {}})
