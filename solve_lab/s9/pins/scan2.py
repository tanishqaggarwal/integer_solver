"""Full pipeline per bit: activate bit + close its pins + re-close mirrors + repair; count failing eqs."""
import sys, time, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/pins')
from build import *
codes,_=H.load_equations()

def act(v, bits):
    s={}
    for b in bits:
        s[b]=1
        for pn in bitpins[b]:
            s[pn['B']]=pn['HUGE']; s[pn['h']]=0
    ripple(v,s)

def pipeline(bits, rounds=6, use8599=True):
    v=list(BASE)
    act(v, bits)
    if use8599 and v[8599]==1 and v[21839]==1:
        ripple(v,{5096:K1, 33612:0})
        ripple(v,{14853:v[12186]})
    close_mirrors(v)
    repair_loop(v, rounds=rounds, verbose=False)
    return v

if __name__=='__main__':
    todo=[b for b in BITS if not BASE[b]]
    if len(sys.argv)>1:
        lo,hi=int(sys.argv[1]),int(sys.argv[2]); todo=todo[lo:hi]
    out={}; t0=time.time()
    for i,b in enumerate(todo):
        v=pipeline([b])
        f=H.evaluate(codes,v)
        out[b]=(len(f), sorted(nz(v)))
        print(f'{b}: fails={len(f)} nz={len(out[b][1])} [{time.time()-t0:.0f}s]',flush=True)
    pickle.dump(out, open(f'pins/scan2_{sys.argv[1] if len(sys.argv)>1 else "all"}.pkl','wb'))
