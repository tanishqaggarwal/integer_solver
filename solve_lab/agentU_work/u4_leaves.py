import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'u_defs.pkl','rb'))
LEAFPIN=D['LEAFPIN']; CONST=D['CONST']
c=collections.Counter(CONST.values())
big=[(v,k) for k,v in c.items() if k.bit_length()>200]
big.sort(reverse=True)
print('big constants in (V-C)/(C-V) pins, top:', [(n, str(k)[:20]+'...', k.bit_length()) for n,k in big[:8]])
# candidate p: the most-repeated 256-bit constant
P=[k for n,k in big if k.bit_length()==256]
print('256-bit constants count', len(P))
for n,k in big[:6]:
    print(n, k.bit_length(), k)
