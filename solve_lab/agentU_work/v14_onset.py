"""V14: what is ON in the deliverable, and is stage 152's coincidence honest or a lie?"""
import pickle, json, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
L=pickle.load(open(B+'v_leaves.pkl','rb')); D=pickle.load(open(B+'v_defs.pkl','rb'))
S=pickle.load(open(B+'v_supp2.pkl','rb')); supp=S['supp']; par=S['par']
T=pickle.load(open(B+'v_tree_final.pkl','rb'))
p=L['p']; shift=L['shift']; sel2exp=L['sel2exp']; pts=L['pts']
V=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
V={int(k[2:]) if k.startswith('x_') else int(k):int(v) for k,v in V.items()}
on=[s for s in sel2exp if V.get(s,0)==1]
print('leaf selectors set to 1 in the deliverable: %d  -> exponents %s'%(len(on),sorted(sel2exp[s] for s in on)))
print('leaf selectors present but not 1:', collections.Counter(V.get(s,0) for s in sel2exp if V.get(s,0)!=1).most_common(3))
# what do the ON leaves' coordinate wires actually carry?
byC={}
for sel,w,C,z,m in D['LEAFPIN']: byC.setdefault(sel,[]).append((w,C%p,m,z))
Xof={(pts[s][0]-shift)%p:sel2exp[s] for s in sel2exp}   # raw u -> exponent
Yof={pts[s][1]:sel2exp[s] for s in sel2exp}
for s in on:
    print(' selector x%d (honest exponent %d):'%(s,sel2exp[s]))
    for w,C,m,z in byC[s]:
        v=V.get(w,0)
        tag = 'HONEST (=own pin C)' if v%p==C else ('carries leaf 2^%s'%Xof.get(v%p, Yof.get(v%p,'?')))
        print('    wire x%-6d m=%-10d value%s  %s   z=x%d -> %s'%(
            w,m,'=0' if v==0 else ' set', tag, z, V.get(z,0)!=0))
# honest fold of the ON set
Nn=115792089237316195423570985008687907852837564279074904382605163141518161494337
tot=sum(1<<sel2exp[s] for s in on)
print('honest scalar of the ON set = %d ; mod N = %d'%(tot, tot%Nn))
# which sibling pair separates the two ON leaves?
E={sel2exp[s] for s in on}
sep=[(I,J) for I,J in T['children'].values() if (E&I) and (E&J)]
print('sibling pairs separating the ON exponents:',[(sorted(I&E),sorted(J&E),len(I),len(J)) for I,J in sep])
