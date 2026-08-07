import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
SEVEN=set([22229,22230,35758,35759,35760,35761,35762])
E12=frozenset([2554,6816,8124,9123,9421,12231,12270,12350,14584,18673,22044,29125])
sub=[a for a in range(L.NA) if frozenset(L.atom2eq[a])<=E12]
print('atoms whose equation set is a SUBSET of E12:', sub)
# equally: atoms adding 0 new equations
print()
# atoms whose eq set equals that of another atom, restricted to the 33 in E12
E12atoms=set()
for e in E12: E12atoms|=set(L.eq_atoms[e][2])
g=collections.defaultdict(list)
for a in sorted(E12atoms): g[frozenset(L.atom2eq[a])].append(a)
print('within E12 atoms, identical-eqset groups of size>1:')
for k,vv in g.items():
    if len(vv)>1: print('   ',vv, 'neq',len(k))
print()
# a19089 vs a19090
for pair in [(19089,19090),(10936,10937),(19087,19088)]:
    a,b=pair
    print(f'eqs(a{a})==eqs(a{b}):', set(L.atom2eq[a])==set(L.atom2eq[b]),
          ' sym-diff', sorted(set(L.atom2eq[a])^set(L.atom2eq[b])))
