import pickle, collections, re, sys
D=pickle.load(open('/home/user/integer_solver/solve_lab/agentU_work/u_atoms.pkl','rb'))
A=D['ATOMS']; EQ=D['EQATOMS']
def shape(n):
    k=n[0]
    if k=='var': return 'V'
    if k=='num': return 'C'
    if k=='neg': return '-'+shape(n[1])
    if k=='add': return '(%s+%s)'%(shape(n[1]),shape(n[2]))
    if k=='sub': return '(%s-%s)'%(shape(n[1]),shape(n[2]))
    if k=='mul': return '(%s*%s)'%(shape(n[1]),shape(n[2]))
c=collections.Counter(shape(v) for v in A.values())
for k,v in c.most_common(40): print('%8d  %s'%(v,k))
print('total shapes',len(c))
# how many equations, how many atoms per equation
print('atoms/eq histogram', collections.Counter(len(e) for e in EQ).most_common(10))
