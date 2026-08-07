"""Support is now identical to the deliverable; find the VALUE difference driving 13 vs 7."""
import sys, json, pickle, collections
from math import gcd
_c2=open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/cansearch2.py').read()
_mark='CODES,_=CK.load_equations()'
_pre=_c2[:_c2.index('print(chr(39)+chr(39))') if False else _c2.index(_mark)]
_pre=_pre[:_pre.rindex(chr(10),0,_pre.rindex('loading checker'))]
_post=_c2[_c2.index('def build2'):_c2.index("if __name__")]
_pre='\n'.join(l for l in _pre.split('\n') if 'import checker as CK' not in l)
exec(_pre)
exec(_post)
vv=build2(24601,2081,27994,vabmode='deliv')
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*NV
for k,x in D.items():
    i=int(k[2:])
    if i<NV: vd[i]=int(x)
r_m=E.run(vv); r_d=E.run(vd)
ATOMS=['((x7075*x8731)+x31864)','((5113045*(x7075*x9118))-x29854)',
       '((x4432-x19964)-x28730)','((x7068-x2099)-(7376877*x642))']
print('atom values (mine vs deliverable):')
for a in ATOMS:
    i=E.residx[a]
    print('  %-52s mine=%-22s deliv=%-22s  ratio_p mine=%s deliv=%s'%(
        a[:52],str(r_m[i])[:20],str(r_d[i])[:20],
        r_m[i]%p==0, r_d[i]%p==0))
print('\nthe four corrupted handle vars:')
for h in (642,28730,29854,31864):
    print('  x%-6d mine=%-24s deliv=%-24s equal=%s'%(h,str(vv[h])[:22],str(vd[h])[:22],vv[h]==vd[h]))
print('\ncofactors u:')
for u in (1329,9413,10903,17325):
    print('  x%-6d mine=%-24s deliv=%-24s equal=%s'%(u,str(vv[u])[:22],str(vd[u])[:22],vv[u]==vd[u]))
# how many wires differ at all now
diff=[w for w in range(NV) if vv[w]!=vd[w]]
print('\ntotal wires differing: %d'%len(diff))
modp=[w for w in diff if (vv[w]-vd[w])%p!=0]
print('  differing mod p      : %d  %s'%(len(modp),modp[:20]))
print('  differing only by a multiple of p: %d'%(len(diff)-len(modp)))
