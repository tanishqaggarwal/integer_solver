"""Re-filter the 378 candidate sites against the baseline-failing equations (M's criterion)."""
import sys, json, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
E=Engine()
M25=[2554,5324,6816,8124,8680,9041,9123,9421,11226,12231,12270,12350,14584,15558,18673,
     21000,22044,22534,22997,28929,29125,29330,32026,35512,38051]
PR=pickle.load(open('price.pkl','rb')); eqs_of={a:set(v) for a,v in PR['eqs_of'].items()}
DELIV=['((x7075*x8731)+x31864)','((5113045*(x7075*x9118))-x29854)',
       '((x4432-x19964)-x28730)','((x7068-x2099)-(7376877*x642))']
# --- CALIBRATION of my incidence map against M's 25 ---
u=set()
for a in DELIV: u|=eqs_of.get(a,set())
print('CALIBRATION: deliverable 4 atoms touch %d equations in my map'%len(u))
print('   rows_target (my map) vs M25 : %d of 25'%len(u&set(M25)))
print('   deliverable checker failures [12231,12270,12350,14584,18673,22044,29125] all in my map? %s'%
      set([12231,12270,12350,14584,18673,22044,29125]).issubset(u))
print('   equations in M25 my map does NOT see:',sorted(set(M25)-u))
rows=json.load(open('candidates.json'))
# reconstruct each row's 4 atoms from its wires
LINK=PR['atom_of']
byvar=collections.defaultdict(list)
from circ2 import vars_of
for a in E.res:
    for v in set(vars_of(E.atoms[a])): byvar[v].append(a)
def atoms_of(r):
    out=[]
    for w in r['slot_wires']:
        if w in LINK: out.append(LINK[w])
    for w in r['vab_wires']:
        c=[a for a in byvar[w] if ('*x%d'%w in a or 'x%d*'%w in a) and a in eqs_of]
        if c: out.append(c[0])
    return out
T=set(M25)
res=[]
for r in rows:
    ats=atoms_of(r)
    uu=set()
    for a in ats: uu|=eqs_of.get(a,set())
    r['rows_target']=len(uu&T); r['atoms']=ats
    res.append(r)
inc=[r for r in res if r['rows_target']>0]
inc.sort(key=lambda r:-r['rows_target'])
print('\nRE-FILTER of all %d candidate sites against the 25:'%len(res))
print('  incident (rows_target > 0): %d'%len(inc))
print('  zero-incident (discard, no pricing needed): %d'%(len(res)-len(inc)))
print('  rows_target histogram:',sorted(collections.Counter(r['rows_target'] for r in res).items(),reverse=True)[:12])
print('\nALL incident sites:')
for r in inc:
    print('  rows_target %-3d handles %-32s c=x%-6d parent x%-6d.%-2s liveleaves %d'%(
        r['rows_target'],','.join('x%d'%h for h in r['handles_h_Pmultiple']),
        r['site_child'],r['parent'],r['side'],r['live_leaves_under']))
json.dump(inc,open('incident_sites.json','w'),indent=0)
# my top-12-by-incidence: confirm M's claim they are 0-incident
top12=sorted(res,key=lambda r:(r['incidence'],-r['depth']))[:12]
print('\nCHECK M\'s claim - my previous top 12 by incidence:')
for r in top12:
    print('   c=x%-6d incidence %-3d  rows_target %d'%(r['site_child'],r['incidence'],r['rows_target']))
