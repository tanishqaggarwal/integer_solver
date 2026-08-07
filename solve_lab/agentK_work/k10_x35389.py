#!/usr/bin/env python3
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
from cascade import Cascade, NV, P
from parse import node_str

C = Cascade()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, val in d.items(): full[int(k[2:])] = int(val)
cand = json.load(open(K + '/cand_k9_handles0.json'))
cv = [0] * NV
for k, val in cand.items(): cv[int(k[2:])] = int(val)
freeset = set(C.E.free)
for w in [35389, 2287, 21889, 25156]:
    print('=== x%d  free=%s  deliverable=%s  cand=%s' % (w, w in freeset, str(full[w])[:70], str(cv[w])[:70]))
    for i in C.var2atoms[w]:
        print('    ', C.atomnames[i][:130], ' | val@cand =', str(C.evala(i, cv))[:50])
