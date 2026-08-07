"""Path shim after the container restart wiped every .pkl in the campaign.

model3.pkl / dag.pkl / orient.pkl were rebuilt from agentE_work/parse3.py and dag.py but
written into agentM_work (I may only write in my own directory).  `harness_m.py` is E's
harness.py with those three paths repointed here.  Importing this module BEFORE `ieng`
registers it under the name `harness`, so E's `import harness as H` resolves to it even
though ieng puts agentE_work first on sys.path.

Nothing about the model changes: parse3 reproduced atoms=40,727 and dag reproduced
free=8,365 / vars-defined=30,383, the same numbers logged before the restart.
"""
import sys, os, importlib.util

MDIR = '/home/user/integer_solver/solve_lab/agentM_work'
os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(20_000_000)

if 'harness' not in sys.modules:
    _spec = importlib.util.spec_from_file_location('harness', os.path.join(MDIR, 'harness_m.py'))
    _h = importlib.util.module_from_spec(_spec)
    sys.modules['harness'] = _h
    _spec.loader.exec_module(_h)

harness = sys.modules['harness']
