#!/usr/bin/env python3
"""Forward engine over Z/mZ for arbitrary modulus m (primes and prime powers) + repair + checkpointing."""
import sys,os,json,pickle,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV,compile_node
from circ2 import vars_of

def egcd(a,b):
    if b==0: return (a,1,0)
    g,x,y=egcd(b,a%b); return (g,y,x-(a//b)*y)
def inv_mod(a,m):
    g,x,_=egcd(a%m,m)
    return None if g!=1 else x%m

class ModM:
    def __init__(self,E,m):
        self.E=E; self.m=m
        src=[]
        for a in E.order:
            c=E.cls[a]
            src.append('v[%d]=(%s)%%m'%(c[1],compile_node(c[2])))
        self.prog=compile('\n'.join(src),'<fm>','exec')
        self.rprog=compile('r[:]=['+','.join('(%s)%%m'%compile_node(E.atoms[a]) for a in E.res)+']','<rm>','exec')
    def run(self,v):
        exec(self.prog,{'v':v,'m':self.m,'__builtins__':{}})
        r=[0]*len(self.E.res)
        exec(self.rprog,{'v':v,'r':r,'m':self.m,'__builtins__':{}})
        return r
    def score(self,r):
        m=self.m
        return [e for e,row in enumerate(self.E.eqres) if sum(k*r[j] for k,j in row)%m]

_ASUP=None
def asup(E):
    global _ASUP
    if _ASUP is None:
        sup=pickle.load(open(os.path.join(HERE,'supp.pkl'),'rb'))
        _ASUP=[]
        for a in E.res:
            s=set()
            for u in vars_of(E.atoms[a]): s|=set(sup[str(u)])
            _ASUP.append(sorted(s))
    return _ASUP

def solve_modm(E,m,bit_a=22106,bit_b=5090,maxit=60,verbose=False):
    M=ModM(E,m); A=asup(E)
    v=[0]*NV; v[bit_a]=1; v[bit_b]=1
    frozen={bit_a,bit_b}
    for it in range(maxit):
        r=M.run(v); nz=[i for i in range(len(r)) if r[i]]
        if verbose: print('  it',it,'nz',len(nz),flush=True)
        if not nz: return v,M,True
        prog=False
        for i in nz:
            r=M.run(v)
            if r[i]==0: continue
            base=r[i]
            for f in A[i]:
                if f in frozen: continue
                old=v[f]; v[f]=(old+1)%m; a1=M.run(v)[i]; v[f]=old
                c=(a1-base)%m
                ic=inv_mod(c,m)
                if ic is None: continue
                v[f]=(old+(-base)*ic)%m
                if M.run(v)[i]==0: frozen.add(f); prog=True; break
                v[f]=old
        if not prog: return v,M,False
    return v,M,False

if __name__=='__main__':
    import sympy
    E=Engine()
    out=os.path.join(HERE,'modm_results')
    os.makedirs(out,exist_ok=True)
    mods=[]
    # 30 primes
    for q in list(sympy.primerange(3,120)): mods.append((q,'prime'))
    for q in [1009,10007,100003,1000003,10000019,1000000007,2**31-1,2**61-1,2**89-1,2**127-1,2**255-19]:
        mods.append((q,'prime'))
    # prime powers
    for q,k in [(2,64),(3,40),(5,30),(7,25),(11,20),(13,20),(1009,8),(65537,4),(1000003,3),(2**31-1,2)]:
        mods.append((q**k,'%d^%d'%(q,k)))
    done=set()
    for m,tag in mods:
        if m in done: continue
        done.add(m)
        f=os.path.join(out,'m_%s.json'%tag.replace('^','p') if tag!='prime' else 'q_%d.json'%m)
        if os.path.exists(f): continue
        t0=time.time()
        try:
            v,M,ok=solve_modm(E,m)
            r=M.run(v); bad=M.score(r)
            nnz=sum(1 for x in r if x)
        except Exception as ex:
            print('m=%s tag=%s ERROR %s'%(m,tag,ex),flush=True); continue
        rec=dict(m=str(m),tag=tag,solved=ok,nonzero_atoms=nnz,failing_eqs=len(bad),secs=round(time.time()-t0,1))
        json.dump(rec,open(f,'w'))
        print('m=%-24s tag=%-8s solved=%-5s nz_atoms=%-4d failing_eqs=%-5d t=%.1f'%(str(m)[:24],tag,ok,nnz,len(bad),time.time()-t0),flush=True)
