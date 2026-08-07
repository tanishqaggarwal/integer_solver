import sys, json, pickle, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
from circ2 import vars_of
E=Engine()
M=pickle.load(open('full_model.pkl','rb')); NODE=M['NODE']; OUT=M['OUT']
rows=json.load(open('incident_atoms.json'))
defrhs={c[1]:c[2] for c in (E.cls[a] for a in E.order)}
# identify the gate variable of each atom and which node it belongs to
sel2node={}
for n in NODE:
    for k in ('sa','sb','sab'):
        for s in NODE[n][k]: sel2node[s]=(n,k)
notof={}
for w,r in defrhs.items():
    if r[0]=='-' and r[1][0]=='c' and r[1][1]==1 and r[2][0]=='v': notof[w]=r[2][1]
for r in rows:
    gates=[]
    for v in set(vars_of(E.atoms[r['atom']])):
        if v in sel2node: gates.append(('sel',v,sel2node[v]))
        elif v in notof and notof[v] in sel2node: gates.append(('NOT',v,sel2node[notof[v]]))
    r['gate']=[(t,'x%d'%v,'x%d'%nd[0],nd[1]) for t,v,nd in gates][:2]
by=collections.defaultdict(list)
for r in rows: by[tuple(sorted(x[2:] for x in r['gate']))].append(r)
print('INCIDENT ATOMS GROUPED BY THE NODE THEY GUARD:')
for k,v in sorted(by.items(),key=lambda kv:-sum(x['rows_target_union'] for x in kv[1])):
    print('  node/gate %-34s  %d atoms, total rt %d'%(str(k),len(v),sum(x['rows_target_union'] for x in v)))
    for r in v: print('       rt %-3d h=x%-6s u=x%-6d %s'%(r['rows_target_union'],r['h'],r['u'],r['atom'][:70]))
DEL={642,28730,29854,31864}
allh=[r['h'] for r in rows]
print('\nCOMPLETE USEFUL SEARCH SPACE: subsets of these %d handles = %d candidates'%(len(allh),2**len(allh)))
print('  handles:',sorted(allh))
print('  the deliverable is the subset',sorted(DEL),'-> 39,026 (calibrated)')
rest=[h for h in allh if h not in DEL]
rt={r['h']:r['rows_target_union'] for r in rows}
cands=[]
for k in range(0,len(rest)+1):
    for c in itertools.combinations(rest,k):
        s=sorted(DEL|set(c)); cands.append((sum(rt[h] for h in s),s))
cands.sort(reverse=True)
print('  supersets of the deliverable set: %d, emitted in descending total rt'%len(cands))
out=dict(incident_atoms=rows,
         all_handles=sorted(allh),
         deliverable_subset=sorted(DEL),
         supersets_of_deliverable=[s for _,s in cands],
         note_baselines=dict(M_baseline_25=True,my_baseline_13=True,
              my13_subset_of_M25=True,filtered_on='union (25)'))
json.dump(out,open('emit_for_M.json','w'),indent=0)
print('\nwrote emit_for_M.json')
for tot,s in cands[:12]: print('   totalrt %-3d %s'%(tot,['x%d'%h for h in s]))
