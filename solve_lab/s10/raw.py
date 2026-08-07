import os, sys, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
v = L.load(os.path.join(LAB,'best','new_instance_partial_39026.json'))
av = L.all_atom_values(v)
f = L.failing_eqs(av)
print('lib failing equation indices :', sorted(f))
print('checker failing line indices :', [12231, 12270, 12350, 14584, 18673, 22044, 29125])
print('MATCH:', sorted(f) == [12231, 12270, 12350, 14584, 18673, 22044, 29125])
print(f'\nNEQ = {L.NEQ}, NA = {L.NA}, NVARS = {L.NVARS}')
sq = sum(1 for i in range(L.NEQ) if L.eq_atoms[i][1])
print(f'squared equations: {sq} of {L.NEQ}')
mult = collections.Counter(L.eq_atoms[i][0] for i in range(L.NEQ))
print(f'multiplier distribution (top): {mult.most_common(8)}')
na = collections.Counter(len(L.eq_atoms[i][2]) for i in range(L.NEQ))
print(f'atoms per equation: {dict(sorted(na.items()))}')
# atoms that are pure pins
pin0 = pin1 = boolean = 0
for a in range(L.NA):
    p = L.polys[a]
    ks = list(p.items())
    if len(ks) == 1 and len(ks[0][0]) == 1 and ks[0][1] in (1, -1): pin0 += 1
    elif len(ks) == 2 and any(len(m) == 2 and m[0] == m[1] for m, c in ks): boolean += 1
    elif len(ks) == 2 and all(len(m) <= 1 for m, c in ks) and any(not m for m, c in ks): pin1 += 1
print(f'\npure "x = 0" atoms: {pin0};  "x = const" atoms: {pin1};  boolean atoms: {boolean}')
