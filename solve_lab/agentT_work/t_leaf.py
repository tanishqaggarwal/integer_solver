#!/usr/bin/env python3
"""AUDIT T30b -- is leaf x34974's residue the SAME shared-wire simultaneity diagnosed at |S|=17,
or a new obstruction?  Load the |S|=8 end state, find the surviving condition, and for every
candidate wire ask: does a root that clears it leave the GLOBAL nonzero count un-decreased
(= blocked by something else it breaks)?  That is exactly closeS4's refusal criterion."""
import os,sys,json,collections,itertools
L='/home/user/integer_solver/solve_lab/agentL_work'; T='/home/user/integer_solver/solve_lab/agentT_work'
os.chdir(L); sys.path.insert(0,L)
g={'__name__':'drv'}
exec(compile(open(os.path.join(T,'from_L','closeS4.py')).read().split("if __name__")[0],'c4','exec'),g)
E=g['E']; SL=g['SL']; SHIFT=g['SHIFT']; p=g['p']; NV=g['NV']
relift=g['relift']; vars_of=g['vars_of']; atomvalvars=g['atomvalvars']
fitc=g['fitc']; roots_c=g['roots_c']; influences=g['influences']; nzcount=g['nzcount']
vv=[0]*NV
for k,val in json.load(open(os.path.join(T,'close_T8.json'))).items(): vv[int(k[2:])]=int(val)
relift(vv); r=E.run(vv)
nz=[E.res[i] for i,x in enumerate(r) if x]
print('|S|=8 end state reloaded: %d nonzero atoms'%len(nz))
for a in nz: print('   %s'%a[:100])
viol=[a for a in SL if r[E.residx[a]]!=0 and SL[a] and r[E.residx[a]]%abs(SL[a])!=0]
print('\nsurviving c>1 violations: %d'%len(viol))
base=nzcount(vv)
print('baseline global nonzero count: %d\n'%base)
for a in viol:
    c=abs(SL[a])//p
    ws=sorted(set(q for q in vars_of(E.atoms[a]) if q in SHIFT)|
              set(q for q in atomvalvars[a] if q in SHIFT))
    print('condition c=%d  atom %s'%(c,a[:64]))
    print('   candidate shift wires: %d'%len(ws))
    blocked=0; noroot=0; helped=0
    for w in ws:
        if not influences(vv,a,w): continue
        C=fitc(vv,a,w,999)
        if C is None: continue
        rs=roots_c(C,c)
        if not rs: noroot+=1; continue
        best=None
        for t in rs[:6]:
            old=vv[w]; vv[w]=old+p*t
            n=nzcount(vv)
            vv[w]=old
            if best is None or n<best: best=n
        if best is not None and best<base: helped+=1; print('      wire x%-6d CLEARS and improves (%d -> %d)'%(w,base,best))
        elif best is not None: blocked+=1; print('      wire x%-6d has a root but global count %d -> %d  BLOCKED'%(w,base,best))
    print('   wires with a root but blocked by collateral: %d ; no root: %d ; helping: %d'%(blocked,noroot,helped))
print('\nINTERPRETATION')
print('  blocked-by-collateral  -> same shared-wire simultaneity as |S|=17: ONE phenomenon.')
print('  no root on every wire  -> a genuinely new obstruction, not simultaneity.')
