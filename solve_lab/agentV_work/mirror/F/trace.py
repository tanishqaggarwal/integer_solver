import sys,pickle
sys.path.insert(0,'.')
from fwd import Engine,NV
from parse import node_str
from circ2 import vars_of
E=Engine()
defmap={}
for a in E.order: defmap[E.cls[a][1]]=a
resby=dict()
for a in E.res:
    c=E.cls[a]
    for u in vars_of(E.atoms[a]): resby.setdefault(u,[]).append(a)
free=set(E.free)
def show(v,depth=0,seen=None):
    if seen is None: seen=set()
    if v in seen or depth>6: return
    seen.add(v)
    pad='  '*depth
    if v in free:
        print(pad+'x%d = FREE  (res atoms: %s)'%(v,[s[:70] for s in resby.get(v,[])]))
        return
    a=defmap[v]
    print(pad+'x%d := %s'%(v,node_str(E.cls[a][2])[:100]))
    for u in sorted(vars_of(E.cls[a][2])): show(u,depth+1,seen)
if __name__=='__main__':
    for t in sys.argv[1:]:
        print('====',t); show(int(t))
