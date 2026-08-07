#!/usr/bin/env python3
import sys,os,pickle,collections,re,time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from parse import parse_line,node_str,const_val
from parse2 import core_of,factors
from circ2 import spine,split_coef,vars_of,evalnode

def classify(a):
    if a[0]=='-':
        if a[1][0]=='v' and a[1][1] not in vars_of(a[2]): return ('def',a[1][1],a[2])
        if a[2][0]=='v' and a[2][1] not in vars_of(a[1]): return ('def',a[2][1],a[1])
    if a[0]=='+':
        if a[1][0]=='v' and a[1][1] not in vars_of(a[2]): return ('def',a[1][1],('neg',a[2]))
        if a[2][0]=='v' and a[2][1] not in vars_of(a[1]): return ('def',a[2][1],('neg',a[1]))
    return ('cons',None,a)

def build():
    t0=time.time()
    lines=open(os.path.join(HERE,'..','..','EQUATIONS.txt')).read().splitlines()
    atoms={}; eqrows=[]; nn=0
    for i,ln in enumerate(lines):
        ln=ln.strip()
        if not ln: continue
        e=parse_line(ln); m,S=core_of(e)
        if S is None: S=e; nn+=1
        row=[]
        for nd in spine(S):
            k,a=split_coef(nd); s=node_str(a)
            if s not in atoms: atoms[s]=a
            row.append((k,s))
        eqrows.append(row)
    print('core_of failed on',nn,'lines'); 
    cls={s:classify(a) for s,a in atoms.items()}
    return atoms,cls,eqrows

if __name__=='__main__':
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
    print('cons shapes:'); 
    for k,v in sh.most_common(15):
        print('  %6d %s'%(v,k[:120]))
    pickle.dump({'atoms':atoms,'cls':cls,'eqrows':eqrows},open(os.path.join(HERE,'circ3.pkl'),'wb'))
