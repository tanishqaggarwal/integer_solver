import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, fast
from gfp import gauss_solve
P = L.P
BITS = fast.BITS
A_CTRL = [16441, 22917, 31339, 33708]
B_CTRL = [13222, 14681, 28486, 38667]
OUTER = [5096, 19750]


def light3(th):
    v = [0] * L.NVARS
    for b in BITS:
        v[b] = 1
    for k, x in th.items():
        v[k] = x
    fw.forward(v)
    for _ in range(3):
        v[14515] = (v[14515] + (v[12186] - v[1308])) % P
        v[21589] = (v[21589] + (v[24908] - v[19083])) % P
        fw.forward(v)
    return v


def blkA(v): return [v[3719] % P, v[25118] % P]
def blkB(v): return [v[25614] % P, v[34220] % P]
def gaps(v): return ((v[12186] - v[1308]) % P, (v[24908] - v[19083]) % P)


def solve_block(th, ctrl, fn, seed, iters=25):
    rnd = random.Random(seed)
    th = dict(th)
    for c in ctrl:
        th[c] = rnd.randrange(1, 1 << 80)
    for it in range(iters):
        v = light3(th)
        r = fn(v)
        if not any(r):
            return th, True
        J = [[0] * len(ctrl) for _ in r]
        for j, c in enumerate(ctrl):
            t2 = dict(th)
            t2[c] = th[c] + 1
            r1 = fn(light3(t2))
            for i in range(len(r)):
                J[i][j] = (r1[i] - r[i]) % P
        d = gauss_solve(J, [(-x) % P for x in r], P)
        if d is None:
            return th, False
        for j, c in enumerate(ctrl):
            th[c] = (th[c] + d[j]) % P
    return th, not any(fn(light3(th)))


if __name__ == '__main__':
    t0 = time.time()
    rnd = random.Random(7)
    th = {c: rnd.randrange(1, 1 << 80) for c in OUTER}
    okA = okB = False
    for seed in range(40):
        if not okA:
            th, okA = solve_block(th, A_CTRL, blkA, seed)
        if not okB:
            th, okB = solve_block(th, B_CTRL, blkB, seed)
        v = light3(th)
        print(f"seed{seed}: A={blkA(v)==[0,0]} B={blkB(v)==[0,0]} gaps={gaps(v)} ({time.time()-t0:.0f}s)", flush=True)
        okA = blkA(v) == [0, 0]
        okB = blkB(v) == [0, 0]
        if okA and okB:
            break
    v = light3(th)
    print("ALL:", blkA(v), blkB(v), gaps(v))
    if blkA(v) == [0, 0] and blkB(v) == [0, 0] and gaps(v) == (0, 0):
        th2 = dict(th)
        th2[14515] = v[14515]
        th2[21589] = v[21589]
        json.dump({str(k): x for k, x in th2.items()}, open('theta_blocks.json', 'w'))
        b = fw.bad_checks(v)
        print(f"SOLVED targets. bad_checks(pre-repair)={len(b)}: {b}")
