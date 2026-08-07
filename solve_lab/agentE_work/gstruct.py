#!/usr/bin/env python3
import pickle, re, ast, collections, json, sys
m=pickle.load(open('model3.pkl','rb'))
atoms=m['atoms']; eqt=m['eq_terms']
VAR=re.compile(r'x_(\d+)')
# classify: atom = "LHS - RHS" at top level with LHS a bare var
defs={}      # var -> list of atom ids defining it
kinds=collections.Counter()
atom_vars=[]
atom_def=[None]*len(atoms)
for i,a in enumerate(atoms):
    vs=set(int(x) for x in VAR.findall(a))
    atom_vars.append(vs)
    t=ast.parse(a,mode='eval').body
    ok=False
    if isinstance(t,ast.BinOp) and isinstance(t.op,ast.Sub) and isinstance(t.left,ast.Name):
        ov=int(t.left.id[2:])
        rv=set(int(x) for x in VAR.findall(ast.unparse(t.right)))
        if ov not in rv:
            atom_def[i]=(ov, ast.unparse(t.right), rv); ok=True
    if ok: kinds['def']+=1
    else: kinds['other']+=1
print("atom kinds", kinds)
for i,d in enumerate(atom_def):
    if d: defs.setdefault(d[0],[]).append(i)
print("vars defined:", len(defs), " multi-def vars:", sum(1 for v in defs.values() if len(v)>1))
mm=collections.Counter(len(v) for v in defs.values())
print("defcount hist", sorted(mm.items()))
# other atoms shapes
oth=collections.Counter()
for i,a in enumerate(atoms):
    if atom_def[i] is None: oth[re.sub(r'x_\d+','V',a)]+=1
print("non-def shapes:", len(oth)); 
for k,v in oth.most_common(15): print("  ",v,repr(k[:120]))
pickle.dump({'atom_def':atom_def,'atom_vars':atom_vars,'defs':defs}, open('struct.pkl','wb'))
