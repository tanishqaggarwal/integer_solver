#!/usr/bin/env python3
"""Forward engine over Z/qZ + a mod-q solver for the whole system."""
import sys,os,json,pickle,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV,compile_node
from circ2 import vars_of

class ModQ:
    def __init__(self,E,q):
        self.E=E; self.q=q
        src=[]
        for a in E.order:
            c=E.cls[a]
            src.append('v[%d]=(%s)%%q'%(c[1],compile_node(c[2])))
        self.prog=compile('\n'.join(src),'<fq>','exec')
        self.rprog=compile('r[:]=['+','.join('(%s)%%q'%compile_node(E.atoms[a]) for a in E.res)+']','<rq>','exec')
    def run(self,v):
        exec(self.prog,{'v':v,'q':self.q,'__builtins__':{}})
        r=[0]*len(self.E.res)
        exec(self.rprog,{'v':v,'r':r,'q':self.q,'__builtins__':{}})
        return r
    def score(self,r):
        q=self.q
        return [e for e,row in enumerate(self.E.eqres) if sum(k*r[j] for k,j in row)%q]

def solve_modq(E,q,bit_a,bit_b,verbose=False):
    M=ModQ(E,q)
    v=[0]*NV; v[bit_a]=1; v[bit_b]=1
    # gauss-seidel repair mod q
    import pickle
    sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
    asup=[]
    for a in E.res:
        s=set()
        for u in vars_of(E.atoms[a]): s|=set(sup[str(u)])
        asup.append(sorted(s))
    frozen={bit_a,bit_b}
    for it in range(60):
        r=M.run(v); nz=[i for i in range(len(r)) if r[i]]
        if verbose: print('  it',it,'nz',len(nz),flush=True)
        if not nz: return v,M,True
        prog=False
        for i in nz:
            r=M.run(v)
            if r[i]==0: continue
            base=r[i]
            for f in asup[i]:
                if f in frozen: continue
                old=v[f]; v[f]=(old+1)%q; a1=M.run(v)[i]; v[f]=old
                c=(a1-base)%q
                if c==0: continue
                v[f]=(old + (-base)*pow(c,q-2,q))%q
                if M.run(v)[i]==0: frozen.add(f); prog=True; break
                v[f]=old
        if not prog: return v,M,False
    return v,M,False

if __name__=='__main__':
    E=Engine()
    for q in [1000003,1000000007,2147483647,2305843009213693951,15485863,104729,65537,999999937,
              4093,7919,1000000000039]:
        t0=time.time()
        v,M,ok=solve_modq(E,q,22106,5090)
        r=M.run(v); bad=M.score(r)
        print('q=%-22d solved=%s  nonzero atoms=%d  equations failing mod q=%d  t=%.1f'%(q,ok,sum(1 for x in r if x),len(bad),time.time()-t0),flush=True)
