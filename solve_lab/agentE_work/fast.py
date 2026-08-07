"""Incremental forward evaluation: only recompute the downstream of changed free vars."""
import sys, collections, math, pickle, os
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
import engine as E
NV=E.NV; SEQ=E.SEQ; definer=E.definer; avars=E.avars; acodes=E.acodes
pos={u:k for k,u in enumerate(SEQ)}
if os.path.exists('users.pkl'):
    users=pickle.load(open('users.pkl','rb'))
else:
    users=collections.defaultdict(list)
    for w in SEQ:
        i,_=definer[w]
        for u in avars[i]:
            if u!=w: users[u].append(w)
    users=dict(users); pickle.dump(users,open('users.pkl','wb'))
atom_of=collections.defaultdict(list)
for i,vs in enumerate(avars):
    for u in vs: atom_of[u].append(i)

def downstream(changed):
    aff=set(); stack=list(changed)
    while stack:
        u=stack.pop()
        for w in users.get(u,()):
            if w not in aff:
                aff.add(w); stack.append(w)
    return aff

def apply_delta(v0, changes):
    """Return (v, affected_vars). v is a fresh copy."""
    v=list(v0)
    for k,val in changes.items(): v[k]=val
    aff=downstream(changes.keys())
    ns={'v':v,'__builtins__':{}}
    for u in sorted(aff, key=lambda u: pos[u]):
        i,kind=definer[u]
        E._solvevar(v,ns,u,i,kind[0])
    return v, aff

def atoms_touching(aff):
    s=set()
    for u in aff: s.update(atom_of[u])
    return s

def resid_delta(v0, base_bad, changes):
    """Return new bad-atom dict, given base_bad for v0."""
    v,aff=apply_delta(v0,changes)
    touched=atoms_touching(set(aff)|set(changes))
    ns={'v':v,'__builtins__':{}}
    bad=dict(base_bad)
    for i in touched:
        r=eval(acodes[i],ns)
        if r: bad[i]=r
        else: bad.pop(i,None)
    return bad, v
