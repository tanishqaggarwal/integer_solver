"""Compare the residual channel partition against the tree's subtree partition at any config."""
import sys,json,collections,random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentM_work')
import mcore as M

T,NODES=M.tree()
ROOT=set(NODES['15298'])
# proper subtree nodes (exclude root), sorted large->small
SUB=[(k,NODES[k]) for k in NODES if k!='15298' and NODES[k]]
SUB.sort(key=lambda kv:-len(kv[1]))

def decompose(S):
    """Greedy maximal-node cover of set S by tree nodes; returns (nodes, leftover)."""
    rem=set(S); used=[]
    for k,g in SUB:
        if g<=rem:
            used.append((k,len(g))); rem-=g
    return used, sorted(rem)

def crossings(part):
    """part: list of (label,set). Return tree nodes split across >1 part."""
    out=[]
    for k,g in SUB:
        hit=[(lab,len(g&s)) for lab,s in part if g&s]
        if len(hit)>1: out.append((k,len(g),hit))
    return out

def analyse(tag, seed, coordfull=True, verbose=True):
    v0,bad0,sig=M.measure(seed,coordfull=coordfull)
    cls=M.classes(sig)
    part=[]
    for k,v in cls.items():
        lab = 'ON' if k=='ON' else 'INERT' if k=='INERT' else 'ch%d'%len(v)
        part.append((lab,set(v)))
    part.sort(key=lambda x:-len(x[1]))
    # unique labels
    seen=collections.Counter()
    p2=[]
    for lab,s in part:
        seen[lab]+=1
        p2.append((lab if seen[lab]==1 else lab+chr(96+seen[lab]), s))
    cr=crossings(p2)
    if verbose:
        print('=== %s   (bad0=%d, classes=%d)'%(tag,len(bad0),len([1 for l,_ in p2 if l.startswith('ch')])))
        for lab,s in p2:
            nd,left=decompose(s)
            print('  %-8s n=%3d  nodes=%s  loose=%d %s'%(lab,len(s),
                  ','.join('%s(%d)'%(a,b) for a,b in nd) or '-',len(left),left if len(left)<=14 else ''))
        print('  CROSSINGS (tree nodes split across classes): %d'%len(cr))
        for k,n,hit in cr[:10]: print('     node %s (%d) -> %s'%(k,n,hit))
    return p2,cr,bad0
