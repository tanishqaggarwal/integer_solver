"""Enumerate the p-quantised solo handles: free inputs occurring in exactly ONE check atom,
   whose effect on that atom is a multiple of p."""
import ev, fast, json
from fast import St, chk
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
FREE=ev.F['free0']
solo=[X for X in FREE if len(chk[X])==1]
print('free inputs occurring in exactly one check atom:',len(solo))
st=St({})
rows=[]; gran={}
for X in solo:
    a=chk[X][0]
    g=st.clone().set_free({X:1})
    d=g.av[a]-st.av[a]
    if d==0: gran['dormant']=gran.get('dormant',0)+1; continue
    if d%p==0:
        gran['p']=gran.get('p',0)+1; rows.append((X,a,d//p))
    else:
        gran['other']=gran.get('other',0)+1
print('granularity census:',gran)
print('p-quantised solo handles:',len(rows))
from collections import Counter
print('distinct carrier atoms:',len(set(a for X,a,m in rows)))
print('sample:',rows[:8])
json.dump([[X,a,str(m)] for X,a,m in rows],open('handles.json','w'))
for nm,X in (('x_7497',7497),('x_11436',11436),('x_22820',22820),('x_14393',14393)):
    hit=[r for r in rows if r[0]==X]
    print('  named handle %s present: %s'%(nm,bool(hit)), hit[:1])
