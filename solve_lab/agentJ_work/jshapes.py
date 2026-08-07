#!/usr/bin/env python3
"""Census of atom shapes."""
import pickle, re, os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
atoms = M['atoms']
eqs = M['eqs']

VAR = re.compile(r'x_\d+')
NUM = re.compile(r'(?<![\w_])\d+')


def shape(a):
    s = VAR.sub('V', a)
    s = NUM.sub('C', s)
    s = s.replace('(', '').replace(')', '')
    return s


c = Counter(shape(a) for a in atoms)
print("distinct shapes:", len(c))
for k, n in c.most_common(60):
    print(f"{n:7d}  {k}")

# how many atoms per equation, and atom reuse
usage = Counter()
for e in eqs:
    for _, aid in e['terms']:
        usage[aid] += 1
print("\natom usage histogram:", Counter(usage.values()).most_common(12))
print("atoms used once:", sum(1 for v in usage.values() if v == 1))
print("total atom slots:", sum(usage.values()))
