"""Find 'pure knob' free variables (free vars whose only atoms are cheap) instance-wide."""
import os, sys, collections, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
P=2**256-2**32-977
S=[22229, 22230, 35758, 35759, 35760, 35761, 35762]
print('atom values in witness:')
for a in S: print(f'  a{a} = {av[a]}   mod p = {av[a]%P}')
free=[u for u in range(L.NVARS) if u not in L.definer]
print('\nfree vars:',len(free))
h=collections.Counter(len(L.var_atoms[u]) for u in free)
print('hist |var_atoms| for free vars:',sorted(h.items())[:12])
# pure knobs: free vars in exactly 1 atom
pure=[u for u in free if len(L.var_atoms[u])==1]
print('free vars in exactly 1 atom:',len(pure))
cnt=collections.Counter()
for u in pure:
    a=L.var_atoms[u][0]; cnt[a]+=1
print('atoms controlled by such knobs:',len(cnt))
fp=collections.Counter(len(L.atom2eq[a]) for a in cnt)
print('footprint hist of knob-controlled atoms:',sorted(fp.items()))
# cluster: group knob-controlled atoms by shared equations -> find dense clusters
KA=sorted(cnt)
print('\nknob-controlled atoms with footprint<=10 :',sum(1 for a in KA if len(L.atom2eq[a])<=10))
json.dump({'pure':pure,'KA':KA},open(os.path.join(HERE,'pa_knobs.json'),'w'))
# which of S are knob-controlled
print('S∩KA =',[a for a in S if a in cnt])
for a in S:
    ks=[u for u in L.avars[a] if u in set(pure)]
    print(f'  a{a} pure knobs: {ks}  freevars: {[u for u in sorted(L.avars[a]) if u not in L.definer]}')
