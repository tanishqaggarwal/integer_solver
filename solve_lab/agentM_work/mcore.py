"""Agent M core: channel measurement at arbitrary configurations, tied to the tree."""
import sys, json, re, collections, itertools, time, pickle, os
E_DIR='/home/user/integer_solver/solve_lab/agentE_work'
F_DIR='/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0,E_DIR)
sys.set_int_max_str_digits(20_000_000)
import engine as E, fast, harness as H
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
ROWS=[7389,10187,20212,20215,28647]
CLUSTERKN=[6083,11436,14393,14853,22820,26489,31339,37012]

def load_seed(name='triple8_seed.json'):
    return {int(k):int(v) for k,v in json.load(open(os.path.join(E_DIR,name))).items()}

def isb(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in ('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)'): return True
    return False

_CAND=None; _BOOLS=None
def bools():
    global _CAND,_BOOLS
    if _BOOLS is None:
        _CAND=sorted(set().union(*[set(E.cone(a)[1]) for a in ROWS]))
        _BOOLS=[f for f in _CAND if isb(f)]
    return _BOOLS

def measure(seed, targets=None, coordfull=False):
    """Flip each boolean 0->1 at this seed; return v0, bad0, {f: signature}."""
    v0=E.forward(seed); bad0=E.badatoms(v0)
    B=targets if targets is not None else bools()
    sig={}
    for f in B:
        if v0[f]!=0:
            sig[f]='ON'; continue
        b1,_=fast.resid_delta(v0,bad0,{f:1})
        d={a:(b1.get(a,0)-bad0.get(a,0)) for a in ROWS}
        if coordfull:
            sig[f]=tuple(x%P for x in (d[7389],d[10187],d[20212],d[20215],d[28647]))
        else:
            sig[f]=((d[20212]+d[28647])%P,(d[20215]+d[10187])%P)
    return v0,bad0,sig

def classes(sig):
    cls=collections.defaultdict(list)
    for f,c in sig.items():
        if c=='ON': cls['ON'].append(f); continue
        if any(c): cls[c].append(f)
        else: cls['INERT'].append(f)
    return cls

# ---- tree ----
def tree():
    t=json.load(open(os.path.join(F_DIR,'tree96.json')))
    return t, {k:set(v['gsup']) for k,v in t.items()}
