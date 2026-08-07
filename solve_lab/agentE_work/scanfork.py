"""Sweep a tree's free bits; each bit runs in a forked child with a hard wall-clock kill,
   so a long C-level linear-algebra call cannot stall the sweep."""
import sys, os, time, pickle, json, signal
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E, bitfeas2 as B
LIM=int(os.environ.get('LIM','45'))
L=pickle.load(open('bitlists.pkl','rb'))
which=sys.argv[1]; bits=L['A'] if which=='A' else L['B']
lo=int(sys.argv[2]) if len(sys.argv)>2 else 0
hi=int(sys.argv[3]) if len(sys.argv)>3 else len(bits)
res={}
for b in bits[lo:hi]:
    rfd,wfd=os.pipe()
    t0=time.time()
    pid=os.fork()
    if pid==0:
        os.close(rfd)
        try:
            base={18956:B.C, b:1}
            bad0,out,S,cols,nonlin,rounds,v0=B.analyse(base,{18956,b})
            maxr,nv,nr,nd,msg,sol=out[-1]
            pins=sorted(set(bad0)-B.CORE); exact=None
            if sol is not None:
                ns=dict(base)
                for f,d in sol.items():
                    if d: ns[f]=v0[f]+d
                v=E.forward(ns); av=E.badatoms(v); ff=E.eqfails(av)
                exact=(len(ff), sorted(set(av)-B.CORE))
                if not exact[1]:
                    json.dump({f"x_{i}":int(v[i]) for i in range(E.NV) if v[i]!=0},
                              open(f'bitsol_{b}_{39033-len(ff)}.json','w'))
            payload=dict(pins=pins,feas=sol is not None,msg=msg,vars=nv,rows=nr,nl=nd,exact=exact,
                         sol={int(k):int(x) for k,x in (sol or {}).items() if x})
        except Exception as e:
            payload=dict(feas=None,msg='ERR '+type(e).__name__)
        os.write(wfd, pickle.dumps(payload)); os.close(wfd); os._exit(0)
    os.close(wfd)
    buf=b''
    import select
    deadline=time.time()+LIM
    while True:
        r,_,_=select.select([rfd],[],[],max(0.1,deadline-time.time()))
        if r:
            c=os.read(rfd,65536)
            if not c: break
            buf+=c
        elif time.time()>deadline:
            os.kill(pid,signal.SIGKILL); break
    os.close(rfd)
    try: os.waitpid(pid,0)
    except Exception: pass
    if buf:
        d=pickle.loads(buf)
        print(f"x_{b}\tpins={d.get('pins')}\tFEAS={d['feas']}\tmsg={str(d.get('msg'))[:70]}\texact={d.get('exact')}\t{time.time()-t0:.0f}s",flush=True)
    else:
        d=dict(feas=None,msg='TIMEOUT')
        print(f"x_{b}\tTIMEOUT\t{time.time()-t0:.0f}s",flush=True)
    res[b]=d
    pickle.dump(res,open(f'scanfork_{which}.pkl','wb'))
