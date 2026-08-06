#!/usr/bin/env python3
"""Find difference-partners of the changed slack vars x_16742, x_14853 in the equations,
so we can co-adjust to heal. Check if partners are free and their equation-multiplicity."""
import json, re
from collections import Counter, defaultdict
from agentA_harness import freeinp, lines, eqvars, NEQ, backward_cone

boolset = set(json.load(open('boolbits.json'))['boolvars'])
# eq multiplicity
allvarcount = Counter()
for i in range(NEQ):
    for x in eqvars[i]: allvarcount[x] += 1

def context(var):
    """Find equations containing x_var and the immediate +/- neighbor pattern."""
    pat = re.compile(r'x_' + str(var) + r'\b')
    diffpat = re.compile(r'\(x_' + str(var) + r'\)\s*-\s*\(x_(\d+)\)')      # (x_var)-(x_y)
    diffpat2 = re.compile(r'\(x_(\d+)\)\s*-\s*\(x_' + str(var) + r'\)')     # (x_y)-(x_var)
    eqs = [i for i in range(NEQ) if var in eqvars[i]]
    partners_after = Counter()   # x_var - x_y
    partners_before = Counter()  # x_y - x_var
    for i in eqs:
        for m in diffpat.finditer(lines[i]): partners_after[int(m.group(1))] += 1
        for m in diffpat2.finditer(lines[i]): partners_before[int(m.group(1))] += 1
    return eqs, partners_after, partners_before

for var, name in [(16742, 'x_16742'), (14853, 'x_14853')]:
    eqs, pa, pb = context(var)
    print(f"=== {name}: appears in {len(eqs)} equations ===")
    print(f"  (x_{var} - x_Y) partners Y: {dict(pa)}")
    print(f"  (x_Y - x_{var}) partners Y: {dict(pb)}")
    # candidate co-adjust partners (appear as x_var - x_Y): to keep diff constant, set x_Y += Δvar
    for y, ct in list(pa.items()) + list(pb.items()):
        isfree = y in freeinp; isbool = y in boolset
        _, stcone = backward_cone(35389); _, ttcone = backward_cone(6671)
        feedsST = (y in stcone) or (y in ttcone)
        print(f"    partner x_{y}: free={isfree} bool={isbool} feedsST={feedsST} in {allvarcount[y]} eqs val-context")

# Also: how does x_16742 appear NON-differically? (e.g. x_27713=C2+x_16742 as sum)
print("\n=== raw occurrences of x_16742 patterns (first few eqs) ===")
cnt = Counter()
for i in range(NEQ):
    if 16742 in eqvars[i]:
        for m in re.finditer(r'.{12}x_16742.{12}', lines[i]):
            cnt[m.group(0)] += 1
for s, c in cnt.most_common(12): print(f"  {c:4d}x  ...{s}...")
