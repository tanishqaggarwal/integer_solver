"""Which of the 15 atoms can EVER be driven by circuit values, and at what selector setting?"""
import sys, json, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentT_work/mirror/F')
from fwd import Engine,NV
E=Engine()
MD=pickle.load(open('full_model.pkl','rb'))
NODE=MD['NODE']; OUT=MD['OUT']; tree=MD['tree']; sub=MD['sub']; ROOT=MD['ROOT']
live=set(MD['live'])
parent={}; side_of={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n; side_of[ch]=s
def lc(n): return len([x for x in sub[n] if x in live])
print('%-8s %-6s %-6s %-6s %-6s  %s'%('node','liveA','liveB','total','selab','verdict'))
for n in (27994,4971,35155,14803,36871):
    a,b=tree[n]
    la,lb=lc(a),lc(b)
    can = la>0 and lb>0
    print('x%-7d %-6d %-6d %-6d %-6s  %s'%(n,la,lb,la+lb,
        'CAN be 1' if can else 'ALWAYS 0',
        'stage checks reachable' if can else 'STAGE CHECKS PERMANENTLY VACUOUS (sel_ab==0 always)'))
print()
print('guard atoms need sel_ab==0 ; stage-check atoms need sel_ab==1 -> at any one node the')
print('two families are MUTUALLY EXCLUSIVE in a single assignment.')
print()
GRP={
 'x27994 vab guards  (need sel_ab=0)':[31864,29854],
 'x27994 stage checks(need sel_ab=1)':[23754,35619,9629],
 'x4971.va slot links(always live)  ':[28730,642],
 'x4971.vb slot link (always live)  ':[37413],
 'x35155 stage checks(need sel_ab=1)':[1844,29305,2892],
 'x14803 stage checks(need sel_ab=1)':[23822,7945],
 'x36871.vb slot links(always live) ':[34113,28355],
}
for k,v in GRP.items(): print('  %-36s %s'%(k,['x%d'%x for x in v]))
DELIV=[642,28730,29854,31864]
print()
print('THE DELIVERABLE = %s'%['x%d'%x for x in DELIV])
print('  -> x29854,x31864 are x27994 GUARDS, so it runs at sel_ab(x27994)=0 (measured: x21279=0).')
print('  -> therefore x23754,x35619,x9629 (that node\'s stage checks) are VACUOUS in it,')
print('     which is why it does not corrupt them - not an oversight.')
