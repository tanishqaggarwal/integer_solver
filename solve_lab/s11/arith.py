import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
P=L.P
for a in [688,1618,40608]:
    print(f"a{a}:")
    for m,c in sorted(L.polys[a].items(), key=lambda kv:(len(kv[0]),kv[0])):
        print("   ", tuple('x%d'%u for u in m), c)
