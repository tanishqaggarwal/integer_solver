"""V1: corrected parse.  Atoms are the MAXIMAL '-' nodes reachable through +, unary -, and *.
Nested differences inside an atom are NOT atoms (that bug merged ua/ub in u1)."""
import sys, re, pickle, collections, time
sys.setrecursionlimit(200000)
SRC='/home/user/integer_solver/EQUATIONS.txt'
B='/home/user/integer_solver/solve_lab/agentU_work/'
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
        while s.peek()=='*': s.i+=1; n=('mul',n,s.factor())
        return n
    def factor(s):
        c=s.peek()
        if c=='-': s.i+=1; return ('neg',s.factor())
        if c=='(':
            s.i+=1; n=s.expr(); assert s.t[s.i]==')'; s.i+=1; return n
        m=re.match(r'x_(\d+)', s.t[s.i:])
        if m: s.i+=m.end(); return ('var',int(m.group(1)))
        m=re.match(r'\d+', s.t[s.i:]); assert m
        s.i+=m.end(); return ('num',int(m.group(0)))
def maxatoms(n,out):
    k=n[0]
    if k=='sub': out.append(n); return
    if k in ('var','num'): return
    if k=='neg': maxatoms(n[1],out); return
    maxatoms(n[1],out); maxatoms(n[2],out)
def s(n):
    k=n[0]
    if k=='var': return 'x%d'%n[1]
    if k=='num': return str(n[1])
    if k=='neg': return '(-%s)'%s(n[1])
    return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(s(n[1]),s(n[2]))
if __name__=='__main__':
    t0=time.time(); EQ=[]; AT={}
    for ln,line in enumerate(open(SRC)):
        line=line.strip()
        if not line: continue
        ast=P(line[:-3].replace(' ','')).expr()
        raw=[]; maxatoms(ast,raw)
        seen=set(); lst=[]
        for a in raw:
            c=s(a)
            if c in seen: continue
            seen.add(c); lst.append(c); AT.setdefault(c,a)
        EQ.append(lst)
    print('equations',len(EQ),'distinct maximal atoms',len(AT),'t=%.1f'%(time.time()-t0))
    print('atoms/eq hist',collections.Counter(len(e) for e in EQ).most_common(6))
    pickle.dump({'EQ':EQ,'AT':AT},open(B+'v_atoms.pkl','wb'))
    def shape(n):
        k=n[0]
        if k=='var': return 'V'
        if k=='num': return 'C'
        if k=='neg': return '-'+shape(n[1])
        return {'add':'(%s+%s)','sub':'(%s-%s)','mul':'(%s*%s)'}[k]%(shape(n[1]),shape(n[2]))
    c=collections.Counter(shape(v) for v in AT.values())
    for k,v in c.most_common(40): print('%8d  %s'%(v,k))
