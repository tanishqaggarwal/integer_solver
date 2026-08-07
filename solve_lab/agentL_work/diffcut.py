"""Diff my generalised cut against the deliverable AT THE DELIVERABLE'S OWN SITE."""
import sys, os, json, pickle, time, collections
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentL_work/cansearch.py').read()
src=src.split("if __name__")[0]
lines=[l for l in src.split(chr(10))
       if 'import checker as CK' not in l and 'CK.load_equations' not in l
       and "loading checker" not in l and 'CK.NVARS' not in l]
src=chr(10).join(lines)
src=src.replace('def exact_fail(vv):','def exact_fail(vv):\n    return None\ndef _unused(vv):')
exec(src)
D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
vd=[0]*NV
for k,x in D.items():
    i=int(k[2:])
    if i<NV: vd[i]=int(x)
rd=E.run(vd)
dset={E.res[i] for i,x in enumerate(rd) if x}
print('deliverable broken atoms (%d):'%len(dset))
for a in sorted(dset): print('   D  ',a[:130])
for sv in (True,False):
    vv=build(24601,2081,27994,set_vab=sv)
    r=E.run(vv)
    mine={E.res[i] for i,x in enumerate(r) if x}
    print('\n--- mine (set_vab=%s): %d broken'%(sv,len(mine)))
    for a in sorted(mine):
        tag='SHARED ' if a in dset else 'EXTRA  '
        print('   %s%s'%(tag,a[:130]))
    print('   deliverable atoms I did NOT break:',[a[:70] for a in sorted(dset-mine)])
    if sv:
        # wire-level diff over every wire my model owns
        own=set([24468,18956])
        for n in NODE:
            for d in OUT[n]: own.update([d['va'],d['vb'],d['vab'],d['out']])
        own.update(M['live']); own.update(M['dead'])
        diff=[w for w in sorted(own) if (vv[w]-vd[w])%p!=0]
        print('\n   wires where I differ from the deliverable mod p: %d'%len(diff))
        for w in diff[:40]:
            role=[]
            for n in NODE:
                for i,d in enumerate(OUT[n]):
                    for k2 in ('va','vb','vab'):
                        if d[k2]==w: role.append('x%d.%s[%d]'%(n,k2,i))
            print('      x%-6d  %-28s mine=%s...  deliv=%s...'%(w,','.join(role[:2]),str(vv[w]%p)[:18],str(vd[w]%p)[:18]))
