#!/usr/bin/env python3
"""AUDIT T29 -- L's closure sweep, run by T.  Scripts handed over in agentT_work/from_L/.
closeS4.py's global guard: accept a shift only if the TOTAL nonzero-atom count strictly
decreases, verified by direct recomputation.
Order: |S|=2 as the CONTROL (must give exactly 2 nonzero atoms of 9,032 and 39,018 via
checker.py -- the number I established in audit T24), then |S| = 3, 5, 8.
Data load needs agentL_work as cwd; the dump must land in agentT_work, so we chdir back.
agentL_work is read-only throughout (PYTHONDONTWRITEBYTECODE=1)."""
import os,sys,json,time,random,traceback
L='/home/user/integer_solver/solve_lab/agentL_work'
T='/home/user/integer_solver/solve_lab/agentT_work'
LAB='/home/user/integer_solver/solve_lab'
os.chdir(L); sys.path.insert(0,L)
g={'__name__':'drv'}
src=open(os.path.join(T,'from_L','closeS4.py')).read().split("if __name__")[0]
t0=time.time()
exec(compile(src,'closeS4','exec'),g)
print('closeS4 prefix loaded (%.0fs)'%(time.time()-t0),flush=True)
close=g['close']; M=g['M']
os.chdir(T)                      # close_<tag>.json now lands in MY directory
sys.path.insert(0,LAB)
import checker as CK
codes,_=CK.load_equations()
print('checker loaded: %d equations'%len(codes),flush=True)
NV=38748
def check(tag):
    fn=os.path.join(T,'close_%s.json'%tag)
    if not os.path.exists(fn): return None,None
    v=[0]*NV
    for k,val in json.load(open(fn)).items(): v[int(k[2:])]=int(val)
    f=CK.evaluate_all(codes,v)
    return len(codes)-len(f), f
PLAN=[('T2ctl',2),('T3',3),('T5',5),('T8',8)]
for tag,n in PLAN:
    rnd=random.Random(7)                       # exactly L's __main__ convention
    S=[24601,2081] if n==2 else rnd.sample(M['live'],n)
    print('\n=== |S|=%d  tag=%s  ON-set=%s ==='%(n,tag,S if n<=8 else '...'),flush=True)
    t=time.time()
    try:
        nz=close(S,tag)
    except Exception:
        print('   *** EXCEPTION ***',flush=True); traceback.print_exc(); continue
    el=time.time()-t
    sc,f=check(tag)
    print('   NONZERO ATOMS = %d of 9032   WALL = %.1fs'%(len(nz),el),flush=True)
    print('   checker: %s  (%s failing)'%(sc,len(f) if f is not None else '?'),flush=True)
    if f is not None and len(f)<=20: print('   failing: %s'%f,flush=True)
    for a in nz[:8]: print('      %s'%a[:104],flush=True)
    if tag=='T2ctl':
        ok = (len(nz)==2 and sc==39018)
        print('   CONTROL %s  (need 2 nonzero atoms and 39018)'%('PASSED' if ok else '*** FAILED ***'),flush=True)
        if not ok:
            print('   stopping per instruction: control did not reproduce the established result',flush=True)
            break
print('\nDONE',flush=True)
