"""WR step 14: from the w=1 branch, take the x_7068 move that clears a22229 and
run the enriched engine on the resulting residual (a29539, a40826)."""
import os, sys
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_frame as W, wr_engine as E
P = ad.P
F = W.F_WIRE
v = L.load(os.path.join(HERE, 'wr_rep1_1_39011.json'))
tgt = T.solve_lin(22229, 7068, v)
v[7068] = tgt
F.fwd(v, rounds=8)
F.report(v, 'after clearing a22229 via x_7068')
T.save(v, os.path.join(HERE, 'wr_w1_x7068.json'))
eng = E.Engine(F, forbid={26064})
eng.run(v, 'w1_x7068', budget=int(sys.argv[1]) if len(sys.argv) > 1 else 2400)
