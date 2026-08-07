#!/usr/bin/env python3
import sys,os,pickle,collections,re,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from parse import parse_line,node_str,const_val
from parse2 import factors
from circ2 import spine,split_coef,vars_of,evalnode
from core2 import core2
from circ3 import classify

def build():
    lines=open(os.path.join(HERE,'..','..','EQUATIONS.txt')).read().splitlines()
    atoms={}; eqrows=[]
    for i,ln in enumerate(lines):
        e=parse_line(ln); m,p,S=core2(e)
        if S is None: S=e
        row=[]
        for nd in spine(S):
            k,a=split_coef(nd); s=node_str(a)
            if s not in atoms: atoms[s]=a
            row.append((k,s))
        eqrows.append(row)
    cls={s:classify(a) for s,a in atoms.items()}
    return atoms,cls,eqrows

if __name__=='__main__':
    t0=time.time()
    atoms,cls,eqrows=build()
    defs=collections.defaultdict(list); cons=[]
    for s,c in cls.items():
        if c[0]=='def': defs[c[1]].append(s)
        else: cons.append(s)
    print('atoms',len(atoms),'defined vars',len(defs),'cons atoms',len(cons))
    print('defs/var hist',sorted(collections.Counter(len(v) for v in defs.values()).items()))
    sh=collections.Counter()
    for s in cons:
        t=re.sub(r'x\d+','X',s); t=re.sub(r'-?\d{2,}','C',t); sh[t]+=1
    for k,v in sh.most_common(20): print('  %6d %s'%(v,k[:110]))
    allv=set()
    for s,a in atoms.items(): allv|=vars_of(a)
    print('vars appearing',len(allv),'max idx',max(allv))
    print('free (never defined) vars',len(allv-set(defs)))
    print('eq len hist',sorted(collections.Counter(len(r) for r in eqrows).items())[:6])
    pickle.dump({'atoms':atoms,'cls':cls,'eqrows':eqrows},open(os.path.join(HERE,'circ4.pkl'),'wb'))
    print('t',time.time()-t0)
