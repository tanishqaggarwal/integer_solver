"""Backward cone tracing over the canonical gate DAG."""
import pickle, collections
import harness as H
NVARS=38748
g = pickle.load(open('gates.pkl','rb'))
polys = pickle.load(open('polys.pkl','rb'))
d = pickle.load(open('atoms.pkl','rb')); src = d['atom_src']
definer, atom_out = g['definer'], g['atom_out']
avars = [set(v for m in P for v in m) for P in polys]
var_atoms = collections.defaultdict(list)
for a,s in enumerate(avars):
    for u in s: var_atoms[u].append(a)

def backcone(roots, maxn=100000):
    seen=set(roots); q=collections.deque(roots); free=set()
    while q:
        u=q.popleft()
        a=definer.get(u)
        if a is None: free.add(u); continue
        for w in avars[a]:
            if w!=u and w not in seen:
                seen.add(w); q.append(w)
                if len(seen)>maxn: return seen, free
    return seen, free

def show(u, depth=0, maxdepth=4, seen=None, v=None):
    if seen is None: seen=set()
    if u in seen or depth>maxdepth: return
    seen.add(u)
    a=definer.get(u)
    pad='  '*depth
    val = '' if v is None else f'  [bits={v[u].bit_length()}]'
    if a is None:
        print(f'{pad}x_{u} = FREE{val}'); return
    print(f'{pad}x_{u} <- [{a}] {src[a][:120]}{val}')
    for w in sorted(avars[a]):
        if w!=u: show(w, depth+1, maxdepth, seen, v)
