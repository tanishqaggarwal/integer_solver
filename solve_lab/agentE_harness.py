#!/usr/bin/env python3
"""Reusable forward-eval harness for the wire-freeing escape.
- Pins all 220 wire members to sign*WIREV (excluded from gate-defining).
- Loads free inputs from a base config (default best_agentD_39018).
- Forward-evaluates all other gates in topo order (re-heals all wiring for any free-input choice).
- Exposes: set_free(v,x), forward(), fails(), L1/L2/L3, and set_core_partners().
Import this; call routines. Running it standalone prints the baseline wire=1 state."""
import json, re, sys
from collections import defaultdict, deque
from agentE_common import build_wire, load_gates, load_lines, p, NVARS, VAR, CORE, NONCORE13
sys.setrecursionlimit(1000000)

M2MOD = 6672769
wire, find2, A = build_wire()
wire_set = set(wire)
gates = load_gates()
lines = load_lines()
eqcode = [compile(VAR.sub(r'v[\1]', L.rsplit('=', 1)[0]), '<e>', 'eval') for L in lines]
eqvars = [set(int(m) for m in VAR.findall(L)) for L in lines]

gate_out = set(t for t, _, _ in gates)
freeinp = set(v for v in range(NVARS) if v not in gate_out)

def build_order(pinned_extra):
    """Topo order of gate targets, excluding wire members and pinned_extra (treated as fixed)."""
    ready = [False]*NVARS
    fixed = wire_set | set(pinned_extra)
    for v in range(NVARS):
        if v not in gate_out or v in freeinp or v in fixed: ready[v] = True
    gu = [0]*len(gates); using = defaultdict(list)
    for gi, (t, rhs, vids) in enumerate(gates):
        u = 0
        for v in vids:
            if not ready[v]: u += 1
            using[v].append(gi)
        gu[gi] = u
    definer = {}; order = []
    q = deque(gi for gi in range(len(gates)) if gu[gi] == 0)
    while q:
        gi = q.popleft(); t, rhs, vids = gates[gi]
        if ready[t] or t in fixed: continue
        definer[t] = gi; order.append(t); ready[t] = True
        for gj in using[t]:
            gu[gj] -= 1
            if gu[gj] == 0: q.append(gj)
    return definer, order

class Harness:
    def __init__(self, base_path='best_agentD_39018.json', wirev=1, pin_partners=True):
        self.wirev = wirev
        base = {int(k[2:]): v for k, v in json.load(open(base_path)).items()}
        self.base = base
        # partners are free inputs we will set explicitly
        self.partners = [30317, 5146, 2936] if pin_partners else []
        self.definer, self.order = build_order(self.partners)
        self.gcode = [compile(VAR.sub(r'v[\1]', gates[self.definer[t]][1]), '<r>', 'eval') for t in self.order]
        self.val = [0]*NVARS
        # init free inputs from base
        for v in freeinp:
            self.val[v] = base.get(v, 0)
        # pin wire
        for w, s in wire.items():
            self.val[w] = s*wirev
        self.ns = {'__builtins__': {}, 'v': self.val}
    def forward(self):
        v = self.val; self.ns['v'] = v
        for k, t in enumerate(self.order):
            v[t] = eval(self.gcode[k], self.ns)
    def set_free(self, var, x):
        self.val[var] = x
    def loads(self):
        return self.val[11150], self.val[25739], self.val[37758]
    def set_core_partners(self):
        L1, L2, L3 = self.loads()
        V = self.wirev
        # x_5101=x_32017=x_26789 = wire sign (all +1) * V  => coefficient V on partner
        # M1 = L1 + x_5101*x_30317 = 0 -> x_30317 = -L1/(sign5101*V)
        s5101 = wire[5101]; s32017 = wire[32017]; s26789 = wire[26789]
        self.val[30317] = -L1 // (s5101*V)
        self.val[2936] = 537773*L3 // (s26789*V)
        assert L2 % (M2MOD) == 0, f"6672769 does not divide L2 (L2%M={L2%M2MOD})"
        # M2 = L2 - 6672769*x_32017*x_5146 = 0 -> x_32017*x_5146 = L2/6672769 -> x_5146 = L2/(6672769*sign*V)
        self.val[5146] = (L2 // M2MOD) // (s32017*V)
        return L1, L2, L3
    def fails(self):
        self.ns['v'] = self.val
        return [i for i in range(len(lines)) if eval(eqcode[i], self.ns) != 0]

if __name__ == '__main__':
    H = Harness(wirev=1)
    H.forward()
    L1, L2, L3 = H.loads()
    print(f"wire=1, base=agentD_39018: L2 % 6672769 = {L2 % M2MOD}")
    F = H.fails()
    core = [i for i in F if i in CORE]; nc = [i for i in F if i not in CORE]
    print(f"before partners: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail); core {len(core)}, noncore {len(nc)}: {sorted(nc)}")
