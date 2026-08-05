import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
p=H.p
# ancestor free-inputs of key gates
def anc(v): return H.anc.get(v,{v})
print("=== ancestry sizes (free inputs feeding each) ===")
for v in [2099,19964,7068,4432,33469,29322,3558,27713,1326,26756,25442,5015,27289,25859]:
    a=anc(v)
    print(f"x_{v}: {len(a)} free ancestors; free? {v in H.freeinp}")
# coupling: shared free ancestors between gaps and core
g_gaps = anc(2099)|anc(19964)|{7068,4432}
g_core = anc(33469)|anc(29322)|anc(3558)|anc(27713)|anc(1326)
print(f"\ngap free-ancestors: {len(g_gaps)}, core free-ancestors: {len(g_core)}")
print(f"shared: {len(g_gaps & g_core)} -> {sorted(g_gaps & g_core)[:20]}")
# the MUX leaves & bits
mux=[6418,12553,9118,8731,31861,14865,2081,4287]
print(f"\nMUX vars in gap-ancestry: {[v for v in mux if v in g_gaps]}")
print(f"MUX vars in core-ancestry: {[v for v in mux if v in g_core]}")
# ancestry of level-1 ripple targets
print("\n=== level-1 ripple ancestry ===")
for v in [26756,25442]:
    a=anc(v)
    inter_mux=[m for m in mux if m in a]
    print(f"x_{v}: {len(a)} free anc, contains 7068/4432? {7068 in a},{4432 in a}, MUX:{inter_mux}")
