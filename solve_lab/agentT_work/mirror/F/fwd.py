#!/usr/bin/env python3
"""Forward evaluator: free inputs -> all vars; then residual constraint / equation scoring."""
import sys,os,pickle,collections,time,json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from circ2 import vars_of
NV=38748

def load():
    d=pickle.load(open(os.path.join(HERE,'circ4.pkl'),'rb'))
    s=pickle.load(open(os.path.join(HERE,'sched.pkl'),'rb'))
    return d,s

def compile_node(n):
    o=n[0]
    if o=='v': return 'v[%d]'%n[1]
    if o=='c': return repr(n[1])
    if o=='neg': return '(-%s)'%compile_node(n[1])
    return '(%s%s%s)'%(compile_node(n[1]),o,compile_node(n[2]))

class Engine:
    def __init__(self):
        d,s=load()
        self.atoms=d['atoms']; self.cls=d['cls']; self.eqrows=d['eqrows']
        self.order=s['order']; self.res=s['usedcons']+s['cons']
        # build assignment program
        src=[]
        for a in self.order:
            c=self.cls[a]
            src.append('v[%d]=%s'%(c[1],compile_node(c[2])))
        self.prog=compile('\n'.join(src),'<fwd>','exec')
        # residual atom program
        self.residx={a:i for i,a in enumerate(self.res)}
        rsrc='r[:]=[' + ','.join(compile_node(self.atoms[a]) for a in self.res) + ']'
        self.rprog=compile(rsrc,'<res>','exec')
        # equation rows in terms of residual atoms only (def atoms are 0)
        self.eqres=[]
        for row in self.eqrows:
            rr=[(k,self.residx[a]) for k,a in row if a in self.residx]
            self.eqres.append(rr)
        self.free=sorted(set(range(NV))-set(c[1] for c in (self.cls[a] for a in self.order)))
    def run(self,v):
        exec(self.prog,{'v':v,'__builtins__':{}})
        r=[0]*len(self.res)
        exec(self.rprog,{'v':v,'r':r,'__builtins__':{}})
        return r
    def score(self,r):
        bad=[]
        for i,rr in enumerate(self.eqres):
            t=0
            for k,j in rr: t+=k*r[j]
            if t: bad.append(i)
        return bad

if __name__=='__main__':
    t0=time.time(); E=Engine(); print('build',time.time()-t0)
    print('free inputs',len(E.free),'residual atoms',len(E.res))
    v=[0]*NV
    t1=time.time(); r=E.run(v); print('fwd',time.time()-t1)
    nz=[i for i,x in enumerate(r) if x]
    print('nonzero residual atoms (all-free=0):',len(nz))
    bad=E.score(r)
    print('failing equations:',len(bad),'=> score',39033-len(bad))
    json.dump({'x_%d'%i:v[i] for i in range(NV) if v[i]},open(os.path.join(HERE,'z0.json'),'w'))
    print('total',time.time()-t0)
