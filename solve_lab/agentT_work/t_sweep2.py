#!/usr/bin/env python3
"""AUDIT T30 -- localise the |S|=8 break.  The sweep's ON-sets are nested:
   T3 = [19745,33287,30242]
   T5 = T3 + [12422,16586]
   T8 = T5 + [35110,3545,34974]
so running |S|=6 and 7 on the SAME prefix pins the break to one added leaf.
Prefix passed EXPLICITLY (not re-drawn from Random(7)) so nesting is guaranteed, not assumed.
   T6 = T5 + [35110]      T7 = T6 + [3545]      (T8 = T7 + [34974])
size horizon  -> the failure appears and worsens with |S|
single leaf   -> exactly one of 35110 / 3545 / 34974 flips a closing run to a failing one."""
import os,sys,json,time,traceback
L='/home/user/integer_solver/solve_lab/agentL_work'
T='/home/user/integer_solver/solve_lab/agentT_work'
LAB='/home/user/integer_solver/solve_lab'
os.chdir(L); sys.path.insert(0,L)
g={'__name__':'drv'}
src=open(os.path.join(T,'from_L','closeS4.py')).read().split("if __name__")[0]
exec(compile(src,'closeS4','exec'),g)
close=g['close']
os.chdir(T)
sys.path.insert(0,LAB)
import checker as CK
codes,_=CK.load_equations()
NV=38748
CHAIN=[19745,33287,30242,12422,16586,35110,3545,34974]
def check(tag):
    fn=os.path.join(T,'close_%s.json'%tag)
    v=[0]*NV
    for k,val in json.load(open(fn)).items(): v[int(k[2:])]=int(val)
    f=CK.evaluate_all(codes,v)
    return len(codes)-len(f), f
print('nested chain: %s'%CHAIN,flush=True)
print('known: |S|=3 CLOSES, |S|=5 CLOSES, |S|=8 FAILS (3 atoms, 39002)',flush=True)
res={}
for n,tag in [(6,'T6'),(7,'T7')]:
    S=CHAIN[:n]
    print('\n=== |S|=%d tag=%s  ON-set=%s   (added leaf vs |S|=%d: %d) ==='%(n,tag,S,n-1,CHAIN[n-1]),flush=True)
    t=time.time()
    try: nz=close(S,tag)
    except Exception:
        print('   *** EXCEPTION ***',flush=True); traceback.print_exc(); continue
    sc,f=check(tag)
    res[n]=(len(nz),sc)
    print('   NONZERO ATOMS = %d of 9032   checker %d (%d failing)   WALL %.1fs'%(len(nz),sc,len(f),time.time()-t),flush=True)
    for a in nz[:8]: print('      %s'%a[:104],flush=True)
print('\n================ VERDICT ================',flush=True)
print('  |S|=3 closes (2 atoms) | |S|=5 closes (2) | |S|=6 %s | |S|=7 %s | |S|=8 FAILS (3, 39002)'%(
      ('%d atoms, %d'%res[6] if 6 in res else '?'),('%d atoms, %d'%res[7] if 7 in res else '?')),flush=True)
if 6 in res and 7 in res:
    c6=res[6][0]==2; c7=res[7][0]==2
    if c6 and c7:   print('  -> SINGLE LEAF: x%d (the 8th) is what breaks it.'%CHAIN[7],flush=True)
    elif c6 and not c7: print('  -> SINGLE LEAF: x%d (the 7th) is what breaks it.'%CHAIN[6],flush=True)
    elif not c6:    print('  -> SINGLE LEAF: x%d (the 6th) is what breaks it.'%CHAIN[5],flush=True)
