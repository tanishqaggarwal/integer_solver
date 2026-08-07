import sys, os, json, time, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2, gs
P=L.P
CAND=[5096,21589,14515,19750,33708,31339,29261,26489,8060,19450,3473,8971,5616,245,
      22162,30213,22820,7497,14393,11436,14853,16742,
      24175,4615,13992,5669,14485,18963,8386,21868,2936,5146,30317]
theta={}
t0=time.time()
v=gs.state(theta)
BAD=[a for a in fw.bad_checks(v)]
print("BAD:", BAD)
v,base,slope=gs.probe(theta, CAND, BAD)
print(f"probe done ({time.time()-t0:.0f}s)")
tab=collections.defaultdict(list)
for (a,c),s in slope.items():
    if s is not None and s%P!=0: tab[a].append(c)
for a in BAD:
    print(f"  a{a}: linear mod-p controls = {tab.get(a,[])}")
json.dump({'BAD':BAD,'tab':{str(k):v_ for k,v_ in tab.items()}}, open('gstab.json','w'))
