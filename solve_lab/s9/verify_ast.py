"""Independent verifier: pure AST walk, no eval/compile/regex-substitution."""
import ast, json, sys, time
NV=38748
def load(path):
    d=json.load(open(path)); v=[0]*NV
    for k,x in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(x)
    return v
def walk(n,v):
    if isinstance(n,ast.Constant): return n.value
    if isinstance(n,ast.Name): return v[int(n.id[2:])]
    if isinstance(n,ast.UnaryOp): return -walk(n.operand,v)
    if isinstance(n,ast.BinOp):
        a=walk(n.left,v); b=walk(n.right,v)
        if isinstance(n.op,ast.Add): return a+b
        if isinstance(n.op,ast.Sub): return a-b
        if isinstance(n.op,ast.Mult): return a*b
    raise ValueError(ast.dump(n))
v=load(sys.argv[1]); t0=time.time(); fails=[]
for i,line in enumerate(open('../EQUATIONS.txt')):
    line=line.strip()
    if not line: continue
    if walk(ast.parse(line.rsplit('=',1)[0],mode='eval').body, v)!=0: fails.append(i)
print(f'[AST verifier] {sys.argv[1]}: {i+1-len(fails)}/{i+1} satisfied, {len(fails)} failing in {time.time()-t0:.0f}s')
print('failing:', fails)
