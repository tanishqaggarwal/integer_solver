"""Verify an assignment with solve_lab/checker.py's own loader+evaluator,
   only raising Python's decimal-string digit cap (values here exceed 4300 digits)."""
import sys, time
sys.set_int_max_str_digits(10_000_000)
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import checker
codes,_=checker.load_equations()
v=checker.load_assignment(sys.argv[1])
t=time.time(); fails=checker.evaluate_all(codes,v)
print(f"[verifyE] {sys.argv[1]}: satisfied {len(codes)-len(fails)}/{len(codes)} ({len(fails)} failing) eval={time.time()-t:.1f}s")
print("failing:",fails[:25])
