"""V11: DIRECT test of the equal-inputs criterion at every gadget, from three gadget-site
sources, with my own supports.  Small gadgets exhaustively; large ones by the exact interval."""
import pickle, json, collections, itertools
B='/home/user/integer_solver/solve_lab/agentU_work/'
S=pickle.load(open(B+'v_supp2.pkl','rb')); supp=S['supp']; par=S['par']
T=pickle.load(open(B+'v_tree_final.pkl','rb'))
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
def find(a):
    while par[a]!=a: a=par[a]
    return a
def mv(s): return sum(1<<e for e in s)
def sp(v): return frozenset(supp.get(find(v),()))

def subsums(s):
    out={0}
    for e in sorted(s):
        w=1<<e; out |= {x+w for x in out}
    return out

def test_pair(I,J,budget=22):
    """returns ('EXH-CLEAR'|'EXH-HIT'|'BOUND-CLEAR'|'OPEN', detail)"""
    a,b=mv(I),mv(J)
    if a<N and b<N: return ('BOUND-CLEAR', max(a,b)/N)
    if len(I)+len(J)<=budget:
        SA=subsums(I); SB=subsums(J)
        for x in SA:
            if (x-N) in SB and (x>0): return ('EXH-HIT',(x,x-N))
            if (x+N) in SB: return ('EXH-HIT',(x,x+N))
        return ('EXH-CLEAR',None)
    return ('OPEN',(a/N,b/N))

def run(name, pairs):
    print('---- %s : %d gadget sites ----'%(name,len(pairs)))
    dis=sum(1 for I,J in pairs if I and J and I.isdisjoint(J))
    print('  slot pairs with disjoint supports: %d/%d'%(dis,len(pairs)))
    res=collections.Counter(); worst=0; hits=[]
    for I,J in pairs:
        st,d=test_pair(I,J); res[st]+=1
        if st=='EXH-HIT': hits.append((I,J,d))
        if st=='BOUND-CLEAR': worst=max(worst,d)
    print('  verdicts:',dict(res),' max maskval/N = %.9f'%worst)
    if hits: print('  *** HITS ***',hits[:3])
    return hits

# ---- source A: my own laminar family (sibling pairs) ----
mine=[(c[0],c[1]) for c in T['children'].values()]
run('U own tree (255 merges)', mine)
# exhaustive brute force on every small sibling pair regardless of the bound
small=[(I,J) for I,J in mine if len(I)+len(J)<=22]
print('  brute-forced exhaustively (|I|+|J|<=22): %d of 255 sibling pairs'%len(small))
cnt=0
for I,J in small:
    SA=subsums(I); SB=subsums(J); cnt+=len(SA)*len(SB)
    assert not any((x-N) in SB or (x+N) in SB for x in SA)
print('  subset-sum pairs enumerated: %d  -- zero representations of +-N'%cnt)

# ---- source B: L's 383-node model, leaf sets mapped through MY exponent map ----
M=pickle.load(open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/full_model.pkl','rb'))
tree=M['tree']; sub=M['sub']
sel2exp=pickle.load(open(B+'v_leaves.pkl','rb'))['sel2exp']
lsup={n:frozenset(sel2exp[l] for l in sub[n] if l in sel2exp) for n in tree}
Lpairs=[(lsup[tree[n][0]],lsup[tree[n][1]]) for n in tree if tree[n] is not None]
run("L's 383-node model (my exponent map)", Lpairs)
print("  L's leaf-set family == my support family :",
      set(lsup.values())-{frozenset()} == set(T['fam']))

# ---- source C: Q's 383 chord gadgets, ua/ub through MY supports ----
st=json.load(open('/home/user/integer_solver/solve_lab/agentQ_work/qstages.json'))['stages']
Qpairs=[(sp(s['ua']),sp(s['ub'])) for s in st]
run("Q's 383 chord gadgets (my supports on ua/ub)", Qpairs)
emp=sum(1 for I,J in Qpairs if not I or not J)
print('  Q sites where a slot has empty support under my map:',emp)
print('  Q distinct unordered pairs:',len({frozenset((I,J)) for I,J in Qpairs}))

# ---- tree-free: EVERY disjoint pair drawn from the 511-set family ----
fam=list(T['fam']); prop=[s for s in fam if len(s)<256]
mx=max(mv(s) for s in prop)
print('---- tree-free ----')
print('  proper supports: %d ; max maskval/N = %.9f ; any >= N : %s'%(len(prop),mx/N,mx>=N))
print('  => for EVERY disjoint pair (I,J) drawn from the family, |sum_A - sum_B| <= max < N,')
print('     so +-N is unreachable without enumerating any pair.')
