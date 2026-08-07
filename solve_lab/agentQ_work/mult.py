import sys, os, json, collections, re
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
sys.path.insert(0, os.path.join(HERE,'..','agentF_work'))
from fwd import Engine, NV
E=Engine()
mult=collections.Counter()
for row in E.eqrows:
    for k,a in row: mult[a]+=1
print('atoms total', len(mult))
print('multiplicity histogram:', sorted(collections.Counter(mult.values()).items())[:15])
print('atoms in exactly 1 equation:', sum(1 for a in mult if mult[a]==1))
PIN1='((x24468-91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002)-x32989)'
PIN2='((8863713*(x18956-125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626))-x14257)'
for nm,a in (('PIN_x',PIN1),('PIN_y',PIN2)):
    print(nm, 'present' if a in mult else 'ABSENT', mult.get(a))
# find them loosely
for a in mult:
    if '91416258160755509' in a or '125787314747601108' in a:
        print('  ->', mult[a], a[:120])
json.dump({a:mult[a] for a in mult}, open('mult.json','w'))
