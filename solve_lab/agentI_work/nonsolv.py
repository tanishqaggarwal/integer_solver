import pickle, os, re, collections
HERE = os.path.dirname(os.path.abspath(__file__))
D = pickle.load(open(os.path.join(HERE,'atoms.pkl'),'rb'))
polys = pickle.load(open(os.path.join(HERE,'polys.pkl'),'rb'))
dag = pickle.load(open(os.path.join(HERE,'dag.pkl'),'rb'))
sol = dag['sol']
groups = collections.defaultdict(list)
for i,s in enumerate(sol):
    if not s:
        t = re.sub(r'X\d+','V', D['atom_src'][i]); t = re.sub(r'\d+','N',t)
        groups[t].append(i)
for t, ids in sorted(groups.items(), key=lambda kv:-len(kv[1])):
    print(f"=== {t}  ({len(ids)})")
    for i in ids[:5]:
        print("     ", D['atom_src'][i], "   poly:", polys[i])
