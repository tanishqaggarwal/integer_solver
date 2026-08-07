"""U10: reconcile my 255-internal-node tree with Q's 383-gadget census and L's 383-node model.
Gadget SITES are taken from Q / L (read-only); the SUPPORTS are mine."""
import pickle, json, collections, sys
B='/home/user/integer_solver/solve_lab/agentU_work/'
S=pickle.load(open(B+'u_supp.pkl','rb')); supp=S['supp']; par=S['par']
def find(a):
    while par[a]!=a: a=par[a]
    return a
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
def mv(s): return sum(1<<e for e in s)
def sp(v):
    return frozenset(supp.get(find(v),()))

print('===== SOURCE 1: Q qstages.json (383 chord gadgets) =====')
st=json.load(open('/home/user/integer_solver/solve_lab/agentQ_work/qstages.json'))['stages']
pairs=[]; empt=0
for s in st:
    I,J=sp(s['ua']), sp(s['ub'])
    if not I or not J: empt+=1
    pairs.append((I,J))
print('383 gadgets; slots with EMPTY support under my map:',empt)
print('distinct unordered (I,J) pairs:', len({frozenset((a,b)) for a,b in pairs}))
print('distinct slot supports used:', len({x for pr in pairs for x in pr}))
dis=sum(1 for a,b in pairs if a and b and a.isdisjoint(b))
print('pairs with disjoint supports: %d/383'%dis)
bad=[(a,b) for a,b in pairs if (a and mv(a)>=N) or (b and mv(b)>=N)]
print('*** gadgets with maskval(slot) >= N :', len(bad))
mx=max(max(mv(a),mv(b)) for a,b in pairs if a and b)
print('max slot maskval / N = %.9f'%(mx/N))
szs=collections.Counter((len(a),len(b)) for a,b in pairs)
print('top (|I|,|J|):', szs.most_common(6))
rootlike=[(a,b) for a,b in pairs if len(a)+len(b)==256]
print('gadgets whose two slots cover all 256:', len(rootlike), set((len(a),len(b)) for a,b in rootlike))

print()
print('===== SOURCE 2: L full_model.pkl (mirrored by T) =====')
try:
    M=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/full_model.pkl','rb'))
    tree=M['tree']; sub=M['sub']
    inter=[n for n in tree if tree[n] is not None]
    print('L nodes:',len(tree),' internal:',len(inter))
    lp=[]
    for n in inter:
        a,b=tree[n]; lp.append((sp(a) if a in supp or find(a) in supp else sp(a), sp(b)))
    # use L's own leaf sets too, mapped through my exponent map
    Ls=pickle.load(open(B+'u_leaves.pkl','rb'))['sel2exp']
    lsup={n:frozenset(Ls[l] for l in sub[n] if l in Ls) for n in tree}
    print('L subtree leaf-sets: distinct =', len(set(lsup.values())))
    prs=[(lsup[tree[n][0]],lsup[tree[n][1]]) for n in inter]
    print('distinct unordered pairs from L:', len({frozenset(x) for x in prs}))
    print('L pairs disjoint: %d/%d'%(sum(1 for a,b in prs if a.isdisjoint(b)), len(prs)))
    badL=[(a,b) for a,b in prs if mv(a)>=N or mv(b)>=N]
    print('*** L gadgets with maskval(slot) >= N :', len(badL))
    print('max L slot maskval / N = %.9f'%(max(max(mv(a),mv(b)) for a,b in prs)/N))
    print('L root split sizes:', sorted((len(a),len(b)) for a,b in prs if len(a)+len(b)==256))
    print('my family == L family ?', set(lsup.values())-{frozenset()} == set(pickle.load(open(B+'u_tree_final.pkl','rb'))['fam']))
except Exception as e:
    print('L model unavailable:', e)
