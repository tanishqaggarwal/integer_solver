"""Clean tree partition: measured from ALL bits off, so no OR-leaf is saturated."""
import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE); sys.path.insert(0, os.path.join(HERE,'..','s9','eff'))
import lib as L
from gmp1 import evalp, forwardp
P=L.P; sys.set_int_max_str_digits(400000)
real=[r[1] for r in json.load(open(os.path.join(HERE,'data','gmp16.json')))]
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp16_base.json')))]
off=list(base)
for b in real: off[b]=0
forwardp(off)
LEAF={8599:'A',21839:'B',7304:'C',25956:'D'}
print("all bits OFF -> leaves:", {n:off[u]%P for u,n in LEAF.items()},
      " U=",off[7715]," V=",off[34554])
grp=collections.defaultdict(list)
for b in real:
    v=list(off); v[b]=1; forwardp(v)
    ch=tuple(n for u,n in LEAF.items() if v[u]!=off[u])
    grp[ch].append(b)
print("\nclean OR-tree partition of the 256 real bits:")
tot=0
for k,vv in sorted(grp.items(), key=lambda z:-len(z[1])):
    print(f"   {k if k else '(drives no leaf)'}: {len(vv)} bits   e.g. {sorted(vv)[:10]}")
    tot+=len(vv)
print("   total",tot)
# U/V side
UV=collections.Counter()
for k,vv in grp.items():
    side = ('U' if set(k)&{'A','B'} else '') + ('V' if set(k)&{'C','D'} else '')
    UV[side or 'none']+=len(vv)
print("\nby MUX side:", dict(UV))
print("\nthe two ON bits:", {b:[n for u,n in LEAF.items() if any(b in vv for kk,vv in grp.items() if n in kk)] for b in (2081,24601)})
json.dump({','.join(k) if k else 'none':sorted(v) for k,v in grp.items()},
          open(os.path.join(HERE,'data','bits_trees.json'),'w'))
