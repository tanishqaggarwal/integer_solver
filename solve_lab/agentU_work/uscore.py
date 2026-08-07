"""Scoring harness: build -> free seed -> M engine forward -> checker.evaluate_all."""
import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
import engine3 as E3, harness as H, checker
ENG=E3.Eng(E3.BASE_DEMOTE)
FREE=set(ENG.FREE)
CODES,_=checker.load_equations()
def seed_of_build(v, extra=None):
    s={k:val for k,val in v.items() if k in FREE and val!=0}
    if extra: s.update({k:val for k,val in extra.items() if val!=0})
    return s
def score(seed):
    vv=ENG.forward(seed)
    return len(checker.evaluate_all(CODES,vv)), vv
