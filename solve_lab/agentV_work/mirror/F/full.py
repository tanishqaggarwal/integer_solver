#!/usr/bin/env python3
"""Incremental full-atom engine over ALL 38,748 variables (no circuit-consistency assumption)."""
import sys,os,pickle,json,collections,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fwd import Engine,NV
from circ2 import vars_of

def compile_node(n):
    o=n[0]
    if o=='v': return 'v[%d]'%n[1]
    if o=='c': return repr(n[1])
    if o=='neg': return '(-%s)'%compile_node(n[1])
    return '(%s%s%s)'%(compile_node(n[1]),o,compile_node(n[2]))

class Full:
    def __init__(self):
        E=Engine()
        self.E=E
        self.names=list(E.atoms)
        self.idx={a:i for i,a in enumerate(self.names)}
        self.NA=len(self.names)
        self.fn=[eval('lambda v: '+compile_node(E.atoms[a])) for a in self.names]
        self.avars=[sorted(vars_of(E.atoms[a])) for a in self.names]
        self.var2atoms=collections.defaultdict(list)
        for i,vs in enumerate(self.avars):
            for u in vs: self.var2atoms[u].append(i)
        self.atom2eq=collections.defaultdict(list)
        self.eqrows=[]
        for e,row in enumerate(E.eqrows):
            rr=[(k,self.idx[a]) for k,a in row]
            self.eqrows.append(rr)
            for k,j in rr: self.atom2eq[j].append((e,k))
        self.NE=len(self.eqrows)
    def init(self,v):
        self.v=v
        self.av=[f(v) for f in self.fn]
        self.ev=[sum(k*self.av[j] for k,j in row) for row in self.eqrows]
        self.nfail=sum(1 for x in self.ev if x)
        return self.nfail
    def setvar(self,i,val):
        """Set v[i]=val, update incrementally, return new nfail."""
        v=self.v; old=v[i]
        if old==val: return self.nfail
        v[i]=val
        nf=self.nfail
        for j in self.var2atoms[i]:
            new=self.fn[j](v); d=new-self.av[j]
            if d==0: continue
            self.av[j]=new
            for e,k in self.atom2eq[j]:
                o=self.ev[e]; n=o+k*d
                if o==0 and n!=0: nf+=1
                elif o!=0 and n==0: nf-=1
                self.ev[e]=n
        self.nfail=nf
        return nf
    def fails(self): return [e for e in range(self.NE) if self.ev[e]]

if __name__=='__main__':
    t0=time.time(); F=Full(); print('build',time.time()-t0)
    d=json.load(open(os.path.join(HERE,'..','best','new_instance_partial_39026.json'))) 
    v=[0]*NV
    for k,x in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(x)
    t1=time.time(); nf=F.init(v); print('init',time.time()-t1,'nfail',nf,'score',39033-nf)
    print('fails',F.fails())
    # timing of a move
    import random
    t2=time.time()
    for _ in range(200):
        i=random.randrange(NV); F.setvar(i,F.v[i]+1); F.setvar(i,F.v[i]-1)
    print('400 moves',time.time()-t2)
