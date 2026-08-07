import sys, os, json, collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
sys.path.insert(0, os.path.join(HERE,'..','agentF_work'))
from fwd import Engine, NV
E=Engine()
DELIV=['((5113045*(x7075*x9118))-x29854)','((x7068-x2099)-(7376877*x642))','((x7075*x8731)+x31864)',
 '(x28730-(x17499*x9413))','(x29854-(x22665*x1329))','(x31864-(x28961*x10903))','(x642-(x28599*x17325))']
FWD=['((x7075*x8731)+x31864)','((5113045*(x7075*x9118))-x29854)','((x4432-x19964)-x28730)','((x7068-x2099)-(7376877*x642))']
for name,SET in (('deliverable-7',DELIV),('forward-4',FWD)):
    S=set(SET); rows=collections.Counter(); forms=collections.Counter()
    for i,row in enumerate(E.eqrows):
        f=tuple(sorted((k,a) for k,a in row if a in S))
        if f: rows[len(f)]+=1; forms[f]+=1
    print('===',name,'atoms',len(S))
    print('  equations touching >=1 of them:', sum(rows.values()), ' by count:', dict(rows))
    print('  distinct restricted forms:', len(forms))
    for f,n in forms.most_common(30):
        print('    x%-3d  %s'%(n, ' + '.join('%d*[%s]'%(k,SET.index(a)) for k,a in f)))
