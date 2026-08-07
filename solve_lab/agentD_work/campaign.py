"""Agent D stochastic campaign.

Single long-running process (checkpointing), designed to survive a kill.
Move set, applied to a state whose gate atoms are re-solved after every move:
  * single-knob exact repair of a nonzero check   (rad knobs, exact division)
  * two-knob Bezout repair                        (closes p-quantised checks)
  * random walk on the handles of currently-nonzero atoms (lattice moves)
  * advice sweep / handle sweep
  * random restarts and perturbations
Acceptance: simulated annealing on -#failing equations, always tracking the best.
Every state with score >= THRESH is written to disk immediately.
"""
import json, sys, time, random, os, math, collections
import dlib as L
import engine2 as E
import rad
import adv3
import hsweep

P = L.P
HERE = os.path.dirname(os.path.abspath(__file__))
THRESH = int(os.environ.get('DTHRESH', '39026'))

occ = collections.Counter()
for a in range(L.NA):
    for u in L.avars[a]:
        occ[u] += 1
SOLO = {u for u in L.freeset if occ[u] == 1}


def knobcache(st, c, cache):
    if c not in cache:
        cache[c] = sorted(rad.free_knobs(c, st.v))
    return cache[c]


def single_repairs(st, c, knobs, rng, limit=80):
    """Exact single-knob solves for atom c."""
    c0 = st.av[c]
    out = []
    ks = list(knobs)
    rng.shuffle(ks)
    for u in ks[:limit]:
        b = st.v[u]
        r = st.apply({u: b + 1})
        s = st.av[c] - c0
        st.revert(r)
        if s == 0 or c0 % s:
            continue
        out.append({u: b - c0 // s})
    return out


def pair_repairs(st, c, knobs, rng, tries=60):
    c0 = st.av[c]
    ks = list(knobs)
    rng.shuffle(ks)
    ks = ks[:24]
    sl = {}
    for u in ks:
        b = st.v[u]
        r = st.apply({u: b + 1})
        s = st.av[c] - c0
        st.revert(r)
        if s:
            sl[u] = s
    items = list(sl.items())
    out = []
    for _ in range(tries):
        if len(items) < 2:
            break
        (u1, s1), (u2, s2) = rng.sample(items, 2)
        g = math.gcd(s1, s2)
        if c0 % g:
            continue
        a0, b0 = _bez(s1, s2)
        k = -c0 // g
        a, b = a0 * k, b0 * k
        m = s2 // g
        if m:
            t = a // m
            a -= t * m
            b += t * (s1 // g)
        out.append({u1: st.v[u1] + a, u2: st.v[u2] + b})
    return out


def _bez(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_s, old_t


def lattice_moves(st, rng, n=14):
    """Random moves on the handles reachable from the currently nonzero atoms:
    shifts the residual inside its coset without changing the placement."""
    nz = st.nz()
    if not nz:
        return []
    out = []
    for _ in range(n):
        c = rng.choice(nz)
        kn = rad.free_knobs(c, st.v)
        cand = [u for u in kn if u in SOLO] or list(kn)
        if not cand:
            continue
        u = rng.choice(cand)
        d = rng.choice([1, -1, 2, -2, 3, -3, 7, -7, rng.randint(-1000, 1000)])
        out.append({u: st.v[u] + d})
    return out


def run(path, seed, budget, tag, useblock=False):
    rng = random.Random(seed)
    v0 = L.load(path)
    blk = ()
    if useblock:
        s_tmp = E.St(v0)
        blk = set(a for a in s_tmp.nz() if a in L.atom_out)
        print(f'[{tag}] blocking gates {sorted(blk)}', flush=True)
    st = E.St(v0, block=blk)
    best = st.score
    bestv = list(st.v)
    print(f'[{tag}] start {path} score={st.score} nz={st.nz()}', flush=True)
    T = 2.0
    t0 = time.time()
    cache = {}
    it = 0
    since = 0
    while time.time() - t0 < budget:
        it += 1
        nz = st.nz()
        if not nz:
            print(f'[{tag}] SOLVED', flush=True)
            L.save(st.v, os.path.join(HERE, 'D_SOLVED.json'))
            return
        moves = []
        c = rng.choice(nz)
        kn = knobcache(st, c, cache)
        moves += single_repairs(st, c, kn, rng)
        moves += pair_repairs(st, c, kn, rng)
        moves += lattice_moves(st, rng)
        if not moves:
            st.apply({rng.choice(sorted(L.freeset)): rng.randint(-3, 3)})
            continue
        scored = []
        for mv in moves:
            r = st.apply(mv)
            scored.append((st.score, mv))
            st.revert(r)
        scored.sort(key=lambda t: -t[0])
        top = scored[:5]
        sc, mv = top[0] if rng.random() > 0.25 else rng.choice(top)
        d = sc - st.score
        if d >= 0 or rng.random() < math.exp(d / max(T, 0.2)):
            st.apply(mv)
            cache.clear()
        since += 1
        if st.score > best:
            best = st.score
            bestv = list(st.v)
            print(f'[{tag}] it{it} NEW BEST {best} nz={st.nz()} t={time.time()-t0:.0f}s', flush=True)
            L.save(bestv, os.path.join(HERE, f'D_camp_{tag}_best.json'))
            if best >= THRESH:
                L.save(bestv, os.path.join(HERE, f'D_CAMP_{best}_{tag}.json'))
        if since > 300:
            since = 0
            T = min(T * 1.5, 20)
            adv3.sweep(st, rounds=4)
            hsweep.sweep(st, rounds=2)
            cache.clear()
        else:
            T = max(T * 0.999, 0.3)
        if it % 200 == 0:
            print(f'[{tag}] it{it} score={st.score} best={best} T={T:.2f} t={time.time()-t0:.0f}s', flush=True)
    print(f'[{tag}] done best={best}', flush=True)
    L.save(bestv, os.path.join(HERE, f'D_camp_{tag}_best.json'))


if __name__ == '__main__':
    path = sys.argv[1]
    seed = int(sys.argv[2])
    budget = float(sys.argv[3])
    tag = sys.argv[4]
    ub = len(sys.argv) > 5 and sys.argv[5] == 'block'
    run(path, seed, budget, tag, ub)
