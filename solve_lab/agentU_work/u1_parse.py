"""U1: independent recursive-descent parse of EQUATIONS.txt -> atom inventory.
No other agent's code or data is read.  Output: u_atoms.pkl
"""
import sys, re, pickle, collections, time
sys.setrecursionlimit(100000)

SRC='/home/user/integer_solver/EQUATIONS.txt'

class P:
    def __init__(s,t): s.t=t; s.i=0
    def peek(s): return s.t[s.i] if s.i<len(s.t) else ''
    def expr(s):
        n=s.term()
        while s.peek() in ('+','-'):
            op=s.t[s.i]; s.i+=1; r=s.term()
            n=('sub',n,r) if op=='-' else ('add',n,r)
        return n
    def term(s):
        n=s.factor()
        while s.peek()=='*':
            s.i+=1; n=('mul',n,s.factor())
        return n
    def factor(s):
        c=s.peek()
        if c=='-':                 # unary minus
            s.i+=1; return ('neg',s.factor())
        if c=='(':
            s.i+=1; n=s.expr()
            assert s.t[s.i]==')', (s.i, s.t[s.i-20:s.i+20])
            s.i+=1; return n
        m=re.match(r'x_(\d+)', s.t[s.i:])
        if m: s.i+=m.end(); return ('var',int(m.group(1)))
        m=re.match(r'\d+', s.t[s.i:])
        assert m, (s.i, s.t[max(0,s.i-40):s.i+40])
        s.i+=m.end(); return ('num',int(m.group(0)))

def atoms_of(n, out):
    """collect every (sub, L, R) node; also recurse."""
    k=n[0]
    if k in ('var','num'): return
    if k=='neg': atoms_of(n[1],out); return
    if k=='sub': out.append(n)
    for c in n[1:]:
        atoms_of(c,out)

def s(n):
    k=n[0]
    if k=='var': return 'x%d'%n[1]
    if k=='num': return str(n[1])
    if k=='neg': return '(-%s)'%s(n[1])
    if k=='add': return '(%s+%s)'%(s(n[1]),s(n[2]))
    if k=='sub': return '(%s-%s)'%(s(n[1]),s(n[2]))
    if k=='mul': return '(%s*%s)'%(s(n[1]),s(n[2]))

def vars_of(n, acc):
    if n[0]=='var': acc.add(n[1]); return
    if n[0]=='num': return
    for c in n[1:]: vars_of(c,acc)

if __name__=='__main__':
    t0=time.time()
    EQATOMS=[]          # per equation: list of canonical atom strings
    ATOMS={}            # canon -> ast
    shapes=collections.Counter()
    for ln,line in enumerate(open(SRC)):
        line=line.strip()
        if not line: continue
        assert line.endswith('= 0'), ln
        body=line[:-3].replace(' ','')
        ast=P(body).expr()
        raw=[]; atoms_of(ast,raw)
        # dedupe by canonical string, keep order
        seen=set(); lst=[]
        for a in raw:
            c=s(a)
            if c in seen: continue
            seen.add(c); lst.append(c); ATOMS.setdefault(c,a)
        EQATOMS.append(lst)
        if ln%5000==0: print(ln, time.time()-t0, flush=True)
    print('equations',len(EQATOMS),'distinct sub-nodes',len(ATOMS),'t=%.1f'%(time.time()-t0))
    pickle.dump({'EQATOMS':EQATOMS,'ATOMS':ATOMS}, open('/home/user/integer_solver/solve_lab/agentU_work/u_atoms.pkl','wb'))
