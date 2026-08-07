"""jm step 6: the exact mod-p Jacobian  d(check)/d(free input)  in frame 2,
for every free input that can reach a pin residue, plus every free input in the
cone of a breakable check.  One forward pass per column.

CHUNKED + RESUMABLE.   usage:  python3 jm_06_jac.py START END
Appends one JSON record per column to jm_jac.jsonl and skips columns already
present, so a killed batch costs at most the columns it had not yet written.
"""
import os, sys, time, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/s10')
import jm_lib as J
import lib as L, tools as T, ad
P = J.P
OUT = '/home/user/integer_solver/solve_lab/s10/jm_jac.jsonl'

PINCONE = [91, 96, 420, 438, 490, 1203, 1530, 1603, 2081, 2441, 2455, 2460, 3203,
           3327, 3420, 3476, 3484, 3545, 3591, 3629, 4012, 4287, 4613, 4701, 4783,
           5090, 5403, 5643, 5910, 6312, 6418, 6467, 6550, 6821, 6882, 7039, 7142,
           7401, 8024, 8370, 8473, 8487, 8731, 8976, 9118, 9260, 9280, 9836, 9859,
           9984, 10422, 10513, 10836, 10852, 11049, 11099, 11189, 11343, 11368,
           12054, 12095, 12553, 13195, 13710, 14273, 14329, 14515, 14808, 14823,
           14853, 14865, 15123, 15880, 15889, 16170, 16348, 16586, 16701, 16827,
           17083, 17378, 17406, 17760, 18022, 18197, 18418, 18623, 18833, 19132,
           19569, 20005, 20007, 20066, 20345, 20414, 20564, 20837, 21074, 21103,
           21112, 21499, 21702, 21743, 21981, 22123, 22351, 22513, 22562, 22623,
           22647, 22730, 22917, 22973, 23262, 23751, 24365, 24367, 24517, 24548,
           25106, 25710, 26005, 26471, 26489, 26907, 27305, 27393, 27711, 28001,
           28005, 28124, 28645, 28704, 28713, 28792, 28955, 29191, 29480, 29785,
           30060, 30448, 30468, 30504, 30505, 31687, 31861, 31899, 32387, 32669,
           32936, 32948, 33177, 33287, 33646, 33708, 34175, 34194, 34246, 34734,
           34955, 34974, 35179, 35740, 35839, 35979, 35981, 36044, 37147, 37236,
           37862, 38055, 38625, 38649]

# checks that any pin move has been seen to break, + the canonical-frame pins
SEEDCHK = [3576, 29539, 40826, 7930, 41512, 21617, 37662, 3578, 19297, 19299,
           30984, 36185, 40812, 3568, 3570, 19088, 22233, 22235, 29465, 36602,
           37887, 688, 1618, 3580, 3582, 25676, 29466, 30976, 30978, 33796,
           40608, 3584, 3586, 7932, 7934, 7936, 29467, 42245, 22231, 31672]


def candidates():
    extra = set(J.DETACH)
    for a in SEEDCHK:
        if a >= L.NA:
            continue
        for u in L.avars[a]:
            extra |= set(x for x in J.cone(u) if x in J.FREESET)
    return sorted(set(PINCONE) | extra)


def make_col(u, VM, CHECKS):
    dv = J.jac_col_full(u, VM)
    out = {}
    for c in CHECKS:
        s = 0
        for x in L.avars[c]:
            dx = dv.get(x, 0)
            if dx:
                s += ad.dpart(c, x, VM) * dx
        if s % P:
            out[c] = s % P
    dR = ((1 if u == 7068 else 0) - dv.get(2099, 0),
          (1 if u == 14853 else 0) - dv.get(1308, 0),
          (1 if u == 24548 else 0) - dv.get(25442, 0))
    dA1 = dv.get(28730, 1 if u == 28730 else 0) % P
    return out, tuple(x % P for x in dR), dA1


def done_set():
    s = set()
    if os.path.exists(OUT):
        for ln in open(OUT):
            try:
                s.add(json.loads(ln)['u'])
            except Exception:
                pass
    return s


if __name__ == '__main__':
    a, b = int(sys.argv[1]), int(sys.argv[2])
    U = candidates()
    print(f'{len(U)} candidates total; this batch [{a},{b})', flush=True)
    W = J.base_state()
    VM = [x % P for x in W]
    CHECKS = J.CHECKS
    have = done_set()
    f = open(OUT, 'a')
    t0 = time.time()
    for i in range(a, min(b, len(U))):
        u = U[i]
        if u in have:
            continue
        c, dR, d1 = make_col(u, VM, CHECKS)
        f.write(json.dumps({'u': u,
                            'col': {str(k): str(v) for k, v in c.items()},
                            'dR': [str(x) for x in dR], 'dA1': str(d1)}) + '\n')
        f.flush()
        if (i - a) % 20 == 0:
            print(f'  {i} x_{u} support {len(c)} ({time.time()-t0:.0f}s)',
                  flush=True)
    f.close()
    print(f'batch done ({time.time()-t0:.0f}s)', flush=True)
