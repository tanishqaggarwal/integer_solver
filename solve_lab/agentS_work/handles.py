"""Census of pure handles over ALL free variables, using structural downstream as a cheap prefilter."""
import sys, collections, pickle, time, json
sys.path.insert(0,'.')
import common as C
import harness as H, engine as E, fast
P=C.P
FREE=sorted(E.FREE)
print("free vars:",len(FREE))
t0=time.time()
struct={}
for f in FREE:
    aff=fast.downstream([f])
    ats=fast.atoms_touching(set(aff)|{f})
    struct[f]=ats
print("structural supports in %.0fs"%(time.time()-t0))
sz=collections.Counter(len(v) for v in struct.values())
print("structural support size histogram (small end):",sorted(sz.items())[:20])
pickle.dump(struct,open('struct_supp.pkl','wb'))
# candidates: structural support of size <= 3
cand=[f for f in FREE if len(struct[f])<=3]
print("candidates |struct sup|<=3:",len(cand))
