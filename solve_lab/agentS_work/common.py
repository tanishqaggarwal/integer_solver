"""Agent S common: load agent E's forward engine (read-only) + helpers."""
import sys, os, json, re, collections, itertools, time, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentE_work')
sys.set_int_max_str_digits(20_000_000)
import harness as H, engine as E, fast, sparse
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
ROWS=[7389,10187,20212,20215,28647]
CLUSTERKN=[6083,11436,14393,14853,22820,26489,31339,37012]
BASE={int(k):int(v) for k,v in json.load(open('triple8_seed.json')).items()}

BOOLPAT=('X - X * X','X * X - X','X * (X - 1)','2 * X * (1 - X)')
def isbool(f):
    for i in H.occ[f]:
        t=re.sub(r'x_%d\b'%f,'X',H.atoms[i])
        if t in BOOLPAT: return True
    return False

def cluster_cone():
    return sorted(set().union(*[set(E.cone(a)[1]) for a in ROWS]))
