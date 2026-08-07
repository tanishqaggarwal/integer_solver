"""(0,1) branch: activate ONE b-tree bit, close the selector atoms with x_14853 / x_18956,
   and solve the whole remaining atom system exactly (no atoms excluded)."""
import sys, time, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B, iterfix
C=B.C
def mux01(s,n=6):
    v=E.forward(s)
    for _ in range(n):
        s[14853]=v[13682]; v=E.forward(v and s)
    return E.forward(s),s
def run(bit, iters=6, log=sys.stdout):
    s={18956:C, bit:1}
    v=E.forward(s)
    for _ in range(6):
        s[14853]=v[13682]; v=E.forward(s)
    ns,hist,ok=iterfix.iterate(s,{18956,bit},iters=iters,exclude=set(),log=log)
    v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
    return len(ff),ns,sorted(av)
if __name__=='__main__':
    for bit in [int(x) for x in sys.argv[1:]]:
        t0=time.time()
        n,ns,av=run(bit)
        print(f"BIT {bit}: fails={n} score={39033-n} bad={av[:10]} ({time.time()-t0:.0f}s)",flush=True)
        v=E.forward(ns)
        json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0}, open(f'f01_{bit}_{39033-n}.json','w'))
        json.dump({str(k):str(int(x)) for k,x in ns.items()}, open(f'f01_{bit}_{39033-n}_seed.json','w'))
