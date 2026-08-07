"""Chunked, checkpointed, resumable sweeps.

Long processes in this sandbox can be killed without warning, so no sweep should
hold results in memory or run for tens of minutes. Pattern:

    from chunk import sweep
    sweep('mysweep', candidates, evaluate, start, end)

- results are appended to s10/runs/<tag>.jsonl the moment each is computed
- a restart skips candidates already recorded, so re-running a batch is free
- drive it from the shell in index ranges:
      for i in 0 200 400 600; do python3 my.py $i $((i+200)); done
"""
import os, sys, json, time

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs')


def done_keys(tag):
    """keys already recorded for this sweep (idempotent resume)."""
    path = os.path.join(RUNS, f'{tag}.jsonl')
    if not os.path.exists(path):
        return set()
    out = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.add(json.loads(line)['key'])
            except Exception:
                continue           # tolerate a truncated final line after a kill
    return out


def sweep(tag, candidates, evaluate, start=0, end=None, keyfn=str, budget=540):
    """Evaluate candidates[start:end], appending each result immediately.

    evaluate(cand) -> dict (or None to skip).  Stops cleanly at `budget` seconds
    so the caller can restart on the next chunk without losing anything.
    """
    os.makedirs(RUNS, exist_ok=True)
    path = os.path.join(RUNS, f'{tag}.jsonl')
    seen = done_keys(tag)
    end = len(candidates) if end is None else min(end, len(candidates))
    t0, n = time.time(), 0
    with open(path, 'a') as f:
        for i in range(start, end):
            c = candidates[i]
            k = keyfn(c)
            if k in seen:
                continue
            if time.time() - t0 > budget:
                print(f'[{tag}] budget reached at index {i}; resume from there',
                      flush=True)
                break
            rec = evaluate(c)
            if rec is None:
                rec = {}
            rec['key'] = k
            rec['i'] = i
            f.write(json.dumps(rec) + '\n')
            f.flush()
            os.fsync(f.fileno())
            n += 1
    print(f'[{tag}] wrote {n} results ({time.time()-t0:.0f}s); '
          f'total recorded {len(done_keys(tag))}', flush=True)
    return path


def load(tag):
    """all results recorded for a sweep so far."""
    path = os.path.join(RUNS, f'{tag}.jsonl')
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except Exception:
                continue
    return out


if __name__ == '__main__':
    # self-test: two chunks, with a simulated kill between them
    cands = list(range(10))
    sweep('selftest', cands, lambda c: {'sq': c * c}, 0, 5)
    sweep('selftest', cands, lambda c: {'sq': c * c}, 5, 10)
    r = load('selftest')
    print(f'self-test: {len(r)} records, resume skipped duplicates: '
          f'{len(r) == len(set(x["key"] for x in r))}')
