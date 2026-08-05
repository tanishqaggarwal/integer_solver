"""Pin catalogue: extract all load pins  G*(x_B - HUGE) - s*x_C  and attach metadata."""
import pickle, re, collections, sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import harness as H
exec(open('repair.py').read().split('if __name__')[0])

P = 2**256 - 2**32 - 977
NV = 38748
d = pickle.load(open('atoms.pkl', 'rb')); src = d['atom_src']
atom2eq = pickle.load(open('atom2eq.pkl', 'rb'))
boolv = set(pickle.load(open('boolvars.pkl', 'rb')))
freeset = set(x for x in range(NV) if x not in definer)

PAT = re.compile(r'^x_(\d+) \* \(x_(\d+) - (\d+)\) - (?:(\d+) \* )?x_(\d+)$')

def catalogue():
    pins = []
    for a, s in enumerate(src):
        m = PAT.match(s)
        if not m: continue
        G, B, HUGE, coef, C = m.groups()
        pins.append(dict(atom=a, G=int(G), B=int(B), HUGE=int(HUGE),
                         coef=int(coef) if coef else 1, C=int(C),
                         neq=len(atom2eq.get(a, [])) if isinstance(atom2eq, dict) else 0))
    return pins

if __name__ == '__main__':
    pins = catalogue()
    print('pin atoms found:', len(pins))
    Gs = collections.Counter(p['G'] for p in pins)
    print('distinct gate vars G:', len(Gs))
    print('G free?:', collections.Counter(('free' if p['G'] in freeset else 'gate') for p in pins))
    print('G boolean?:', collections.Counter((p['G'] in boolv) for p in pins))
    print('B free?:', collections.Counter(('free' if p['B'] in freeset else 'gate') for p in pins))
    print('C free?:', collections.Counter(('free' if p['C'] in freeset else 'gate') for p in pins))
    print('distinct B:', len(set(p['B'] for p in pins)), ' distinct C:', len(set(p['C'] for p in pins)))
    print('top G:', Gs.most_common(12))
    print('coef histogram(top):', collections.Counter(p['coef'] for p in pins).most_common(8))
    pickle.dump(pins, open('pins/pins.pkl', 'wb'))
