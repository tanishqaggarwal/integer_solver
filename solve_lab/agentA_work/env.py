"""Agent A shared environment: load lib, frame2 forward, timings."""
import os, sys, time, json
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, os.path.join(LAB, 's10'))
import lib as L
import tools as T
P = 2**256 - 2**32 - 977
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')

def load_best():
    return L.load(BEST)
