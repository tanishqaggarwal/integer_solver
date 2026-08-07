#!/usr/bin/env python3
"""AUDIT T2: agent I's eq8680 hunt tested 22 of the 27 nonzero-net-effect knob groups
in its OWN census (eq8680.log).  Run the 5 it never tested."""
import os,sys,time,json,itertools
I=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','agentI_work'))
sys.path.insert(0,I); os.chdir(I)
from eq8680 import build, minfail_bnb, M, wit, av, v2a, BASE, RES
E_BASE=sorted({e for a in BASE for e,_ in M.atom_eqs[a]}); SB=set(E_BASE)
MISSING=[('X21279',(1629,19227,20291,23436,23438,35829,36502,36505)),
         ('X23754',(23435,23436)),
         ('X34600',(11773,18074,22330,34813,38498)),
         ('X3629',(11921,13644,20293,30859)),
         ('X8976',(11915,11927,20295))]
def show(tag,add,tl):
    SUP=BASE+[a for a in add if a not in BASE]
    r=build(SUP)
    if r is None:
        print('  %s add=%s: non-linear knob -> SKIPPED by build()'%(tag,add),flush=True); return
    SUPs,knobs,E,base,Mat=r
    new=[e for e in E if e not in SB]; t0=time.time()
    mf,forced,nact,nodes=minfail_bnb(E,base,Mat,budget=6,tlimit=tl)
    txt='minfail > 6  (cannot beat 39,026)' if mf is None else ('TIMEOUT %d nodes'%nodes if mf=='timeout' else 'minfail = %s'%mf)
    star='   *** BEATS 39,026 ***' if isinstance(mf,int) and mf<7 else ''
    print('  %-8s add=%-46s |E|=%3d (+%2d new) knobs=%2d forced=%s active=%s nodes=%6d %4ds  %s%s'
          %(tag,str(add),len(E),len(new),len(knobs),forced,nact,nodes,int(time.time()-t0),txt,star),flush=True)
tl=int(sys.argv[1]) if len(sys.argv)>1 else 900
print('base support |E|=%d'%len(E_BASE),flush=True)
print('--- 5 nonzero-effect groups agent I never tested ---',flush=True)
for tag,add in MISSING: show(tag,add,tl)
print('--- pairs of those 5 with the cheapest compensator X19964 (net -1) ---',flush=True)
for tag,add in MISSING: show('X19964+'+tag,tuple(sorted(set(add)|{1631})),tl)
print('done',flush=True)
