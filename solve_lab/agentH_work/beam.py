"""Beam search over pin-cascade repairs in the agent-H frame."""
import ev, engine as E, json, os, sys, time, random
HERE=os.path.dirname(os.path.abspath(__file__))

def key(fv):
    v,f,nz=E.state(fv)
    return (len(f),len(nz)),v,f,nz

def moves(fv,nz):
    out=[]
    for a in nz:
        for u in E.SUPV[a]:
            r=E.solve_u(fv,a,u)
            if r=='ANY' or not r: continue
            for rt in r:
                if fv.get(u,0)==rt: continue
                out.append((a,u,rt))
    return out

def beam(fv0, width=24, depth=60, log=None, seed=0):
    rnd=random.Random(seed)
    k0,v,f,nz=key(fv0)
    frontier=[(k0,fv0)]
    seen=set()
    best=(k0,dict(fv0))
    for d in range(depth):
        cand=[]
        for k,fv in frontier:
            _,_,f_,nz_=key(fv)
            mv=moves(fv,nz_)
            for (a,u,rt) in mv:
                g=dict(fv); g[u]=rt
                sig=tuple(sorted(g.items()))
                if sig in seen: continue
                seen.add(sig)
                kk,_,ff,nn=key(g)
                cand.append((kk,g))
        if not cand: break
        cand.sort(key=lambda t:t[0])
        frontier=cand[:width]
        if frontier[0][0]<best[0]:
            best=(frontier[0][0],dict(frontier[0][1]))
        msg='depth %d: best (fail,nzatoms)=%s  score=%d  cands=%d'%(d,frontier[0][0],39033-frontier[0][0][0],len(cand))
        print(msg,flush=True)
        if log: log.write(msg+'\n'); log.flush()
        if best[0][0]==0: break
    return best

if __name__=='__main__':
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    ub=int(sys.argv[1]) if len(sys.argv)>1 else 542
    wb=int(sys.argv[2]) if len(sys.argv)>2 else 1438
    fv={5096:C2,21589:C1,16742:C1,12186:C2,18956:C1,24468:C2,ub:1,wb:1}
    t0=time.time()
    log=open(os.path.join(HERE,'beam_%d_%d.log'%(ub,wb)),'w')
    (k,bfv)=beam(fv,width=int(os.environ.get('W','24')),depth=int(os.environ.get('D','60')),log=log)
    print('BEST',k,'score',39033-k[0],'time %.1f'%(time.time()-t0))
    v=ev.forward(bfv)
    out=os.path.join(HERE,'beam_%d_%d.json'%(ub,wb))
    json.dump({('x_%d'%i):v[i] for i in range(38748) if v[i]!=0},open(out,'w'))
    json.dump({str(a):b for a,b in bfv.items()},open(out.replace('.json','_fv.json'),'w'))
