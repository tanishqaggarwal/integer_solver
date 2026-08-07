#!/usr/bin/env python3
"""Text-level (paren-accurate) parse of EQUATIONS.txt into atoms. Agent J.

Grammar observed:
  LINE  := LHS ' = 0'
  LHS   := (C)*(CORE) | (C1)*(CORE)+(C2)*(CORE) | (CORE)*(CORE) | (CORE)
           | (C)*((-1)*(CORE))
  CORE  := ((CORE)+(c)*(ATOM)) | (ATOM0)
The chain suffix is literally  '+(<int>)*(' ... ')'  so we peel from the right
by paren matching, requiring the coefficient group to be a bare integer.
"""
import re, sys, pickle, time, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
INT_RE = re.compile(r'^-?\d+$')


def match_fwd(s, i):
    """s[i]=='(' -> index of matching ')'."""
    d = 0
    for j in range(i, len(s)):
        c = s[j]
        if c == '(':
            d += 1
        elif c == ')':
            d -= 1
            if d == 0:
                return j
    raise ValueError('unbalanced')


def strip_parens(s):
    s = s.strip()
    while s.startswith('(') and match_fwd(s, 0) == len(s) - 1:
        s = s[1:-1].strip()
    return s


def peel_chain(s):
    """CORE text -> list of (coef, atom_text), left to right."""
    terms = []
    cur = strip_parens(s)
    while True:
        if not cur.startswith('('):
            break
        m = match_fwd(cur, 0)
        rest = cur[m + 1:]
        if not rest or rest[0] not in '+-':
            break
        sgn = 1 if rest[0] == '+' else -1
        r = rest[1:]
        if not r.startswith('('):
            break
        m2 = match_fwd(r, 0)
        cf = r[1:m2].strip()
        if not INT_RE.match(cf):
            break
        r2 = r[m2 + 1:]
        if not r2.startswith('*('):
            break
        if match_fwd(r2, 1) != len(r2) - 1:
            break
        atom = r2[2:-1]
        terms.append((sgn * int(cf), atom))
        cur = strip_parens(cur[1:m])
    # base
    terms.append((1, cur))
    terms.reverse()
    return terms


def split_top(s, ops='+-'):
    """split at top-level (depth 0) + / - ; returns list of (sign, piece)."""
    out = []
    d = 0
    start = 0
    sign = 1
    for j, c in enumerate(s):
        if c == '(':
            d += 1
        elif c == ')':
            d -= 1
        elif d == 0 and c in ops and j > 0:
            out.append((sign, s[start:j]))
            sign = 1 if c == '+' else -1
            start = j + 1
    out.append((sign, s[start:]))
    return out


def as_int_mult(s):
    """'(int)*(REST)' -> (int, REST) at top level, else None."""
    s = s.strip()
    if not s.startswith('('):
        return None
    m = match_fwd(s, 0)
    head = s[1:m].strip()
    rest = s[m + 1:]
    if rest.startswith('*(') and INT_RE.match(head) and match_fwd(rest, 1) == len(rest) - 1:
        return int(head), rest[1:]
    return None


def outer(s, depth=0):
    """LHS -> (kind, mult, core_text)."""
    s = strip_parens(s)
    parts = split_top(s)
    # (C1)*(S)+(C2)*(S)
    if len(parts) == 2:
        a = as_int_mult(parts[0][1])
        b = as_int_mult(parts[1][1])
        if a and b:
            ca = strip_parens(a[1]); cb = strip_parens(b[1])
            if ca == cb:
                return 'lin', parts[0][0] * a[0] + parts[1][0] * b[0], ca
        return 'lin', 1, s
    # (C)*(REST)
    im = as_int_mult(s)
    if im is not None:
        k, m, c = outer(im[1], depth + 1)
        return k, m * im[0], c
    # (S)*(S)
    if s.startswith('('):
        m = match_fwd(s, 0)
        head = s[1:m]
        rest = s[m + 1:]
        if rest.startswith('*(') and match_fwd(rest, 1) == len(rest) - 1:
            th = strip_parens(rest[1:])
            hh = strip_parens(head)
            if hh == th:
                return 'sq', 1, hh
    return 'lin', 1, s


def main():
    t0 = time.time()
    eqs = []
    atom_ids = {}
    atoms = []
    shapes = Counter()
    with open(EQ) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            lhs = line.rsplit('=', 1)[0].strip()
            kind, mult, core = outer(lhs)
            lst = []
            for c, a in peel_chain(core):
                a = strip_parens(a)
                aid = atom_ids.get(a)
                if aid is None:
                    aid = len(atoms); atom_ids[a] = aid; atoms.append(a)
                lst.append((c, aid))
            eqs.append({'i': idx, 'kind': kind, 'mult': mult, 'terms': lst})
            shapes[kind] += 1
            if idx % 10000 == 0:
                print(f"  {idx} ... {time.time()-t0:.1f}s", file=sys.stderr)
    print(f"parsed {len(eqs)} eqs, {len(atoms)} distinct atoms in {time.time()-t0:.1f}s")
    print("kinds:", dict(shapes))
    with open(os.path.join(HERE, 'jmodel2.pkl'), 'wb') as f:
        pickle.dump({'eqs': eqs, 'atoms': atoms}, f)
    for a in atoms[:25]:
        print("ATOM:", a)


if __name__ == '__main__':
    main()
