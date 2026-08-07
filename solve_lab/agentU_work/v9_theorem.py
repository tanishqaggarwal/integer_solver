"""U9: the partition theorem, tested directly on the recovered support family."""
import pickle, collections, itertools
B='/home/user/integer_solver/solve_lab/agentU_work/'
S=pickle.load(open(B+'v_supp.pkl','rb'))
fam=list(S['fam'])
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
def mv(s): return sum(1<<e for e in s)
fam.sort(key=len)
print('family size',len(fam))
# ---- laminarity + binary-tree reconstruction ----
byset={s:i for i,s in enumerate(fam)}
bad=0
for a,b in itertools.combinations(fam,2):
    if a&b and not (a<=b or b<=a): bad+=1
print('laminarity violations:',bad)
children={}
for s in fam:
    if len(s)==1: continue
    subs=[t for t in fam if t<s]
    maximal=[t for t in subs if not any(t<u<s for u in subs)]
    children[s]=maximal
sz=collections.Counter(len(v) for v in children.values())
print('children-count histogram over internal nodes:', sorted(sz.items()))
print('all internal nodes binary & partition:',
      all(len(v)==2 and v[0].isdisjoint(v[1]) and v[0]|v[1]==k for k,v in children.items()))
root=max(fam,key=len)
print('root size',len(root),'== {0..255}:', root==frozenset(range(256)))
A,Bh=children[root]
if len(A)<len(Bh): A,Bh=Bh,A
print('root halves: %d / %d  disjoint=%s'%(len(A),len(Bh),A.isdisjoint(Bh)))
hi=set(range(129,256))
print('  A omits %d exponents >=129 ; B omits %d'%(len(hi-A),len(hi-Bh)))
print('  |A & hi| = %d, |B & hi| = %d'%(len(A&hi),len(Bh&hi)))
# ---- THE TEST: maskval >= N is necessary for the +N side ----
slots=[s for s in fam if s!=root]
print('proper slot supports (candidate sides):',len(slots))
over=[s for s in slots if mv(s)>=N]
print('supports with maskval >= N :', len(over))
print('max maskval over proper supports / N = %.6f'%(max(mv(s) for s in slots)/N))
top=sorted(slots,key=lambda s:-mv(s))[:5]
for s in top:
    print('   size %3d  maskval/N = %.9f  contains {129..255}? %s  missing>=129: %d'%(
        len(s), mv(s)/N, hi<=s, len(hi-s)))
# ---- and the both-sides check per actual gadget ----
viol=[]
for k,(c1,c2) in children.items():
    if mv(c1)>=N or mv(c2)>=N: viol.append(k)
print('gadgets where either slot could reach +-N (maskval>=N):', len(viol))
# ---- adjacent hole: a half folding to the identity needs sum_S 2^e == N exactly ----
print('supports S with maskval(S) >= N (needed for sum_S = N):', len(over))
suppN=set(e for e in range(256) if (N>>e)&1)
print('any support containing supp(N)?', any(suppN<=s for s in slots))
pickle.dump({'fam':fam,'children':children,'root':root}, open(B+'v_tree_final.pkl','wb'))
