"""Reverse-engineer circuit SEMANTICS: print the definition DAG upward from a var."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L

P = L.P
def fmt_poly(a):
    Pp = L.polys[a]
    parts=[]
    for m,c in sorted(Pp.items(), key=lambda kv:(len(kv[0]), kv[0])):
        s = ('%+d'%c) if c not in (1,-1) or not m else ('+' if c==1 else '-')
        if m: s += '*'.join('x%d'%u for u in m)
        else: s = '%+d'%c
        parts.append(s)
    return ' '.join(parts)

def const_name(c):
    if c==P: return 'P'
    if c==-P: return '-P'
    return str(c)

def show(var, depth, seen, indent=0, maxdepth=6):
    pad='  '*indent
    d = L.definer.get(var)
    if d is None:
        print(f"{pad}x{var} = FREE")
        return
    print(f"{pad}x{var} := a{d}: {fmt_poly(d)}")
    if indent>=maxdepth: return
    for u in sorted(L.avars[d]):
        if u==var: continue
        if u in seen: 
            continue
        seen.add(u)
        show(u, depth, seen, indent+1, maxdepth)

if __name__=='__main__':
    for v in [int(x) for x in sys.argv[1:]]:
        print("="*70); show(v, 0, set(), 0, int(os.environ.get('D','4')))
