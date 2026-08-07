#!/usr/bin/env python3
"""Q-13: does satisfying the bundled equations FORCE the individual atoms to zero?
Each equation is a weighted sum of atoms.  If an atom occurs in only one equation, that equation
alone cannot force it to vanish.  Count, WITHOUT deduplication, how many distinct equations each
atom term occurs in."""
import ast,re,collections,json,pickle
def const_val(n):
    if isinstance(n,ast.Constant): return n.value
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub) and isinstance(n.operand,ast.Constant): return -n.operand.value
    return None
def strip_outer(n):
    while True:
        if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Mult):
            a,b=n.left,n.right; ca,cb=const_val(a),const_val(b)
            if ca is not None and cb is not None: return n
            if ca is not None: n=b; continue
            if cb is not None: n=a; continue
            if ast.unparse(a)==ast.unparse(b): n=a; continue
            return n
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): n=n.operand; continue
        return n
def flatten(n):
    out=[]
    def rec(x):
        if isinstance(x,ast.BinOp) and isinstance(x.op,ast.Add): rec(x.left); rec(x.right)
        else: out.append(x)
    rec(n); return out
occ=collections.defaultdict(set); per_eq=[]
for ei,L in enumerate(open('/home/user/integer_solver/EQUATIONS.txt')):
    L=L.strip()
    if not L: continue
    node=ast.parse(L.rsplit('=',1)[0].strip(),mode='eval').body
    ts=[ast.unparse(strip_outer(t)) for t in flatten(strip_outer(node))]
    per_eq.append(len(ts))
    for s in ts: occ[s].add(ei)
mult=collections.Counter(len(v) for v in occ.values())
print('distinct atom terms: %d'%len(occ))
print('atoms by number of DISTINCT equations containing them:')
for k in sorted(mult): print('   %3d equations : %6d atoms'%(k,mult[k]))
tot=sum(mult.values()); ge2=sum(v for k,v in mult.items() if k>=2)
print('atoms occurring in >= 2 equations: %d / %d  (%.1f%%)'%(ge2,tot,100.0*ge2/tot))
import statistics
print('atoms per equation: mean %.1f  min %d  max %d'%(statistics.mean(per_eq),min(per_eq),max(per_eq)))
ST=[x for x in json.load(open('qstages.json'))['stages'] if 'u3' in x]
core=set()
for s in ST: core|={s['R1'],s['S'],s['dx'],s['dy']}
VAR=re.compile(r'^x_(\d+) -')
cm=collections.Counter()
for s,eqs in occ.items():
    m=VAR.match(s)
    if m and int(m.group(1)) in core: cm[len(eqs)]+=1
print('stage-core atoms by equation multiplicity:',dict(sorted(cm.items())))
pickle.dump({k:len(v) for k,v in occ.items()},open('qmult.pkl','wb'))
