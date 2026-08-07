"""Agent B: classify the 39k gate factors and recover the circuit."""
import pickle, collections, random, sys

W = '/home/user/integer_solver/solve_lab/agentB_work/'
M = pickle.load(open(W+'model5.pkl','rb'))
facs, atoms, eqs = M['facs'], M['atoms'], M['eqs']

def shape(p):
    """canonical shape string, with variables anonymised"""
    ms = sorted(p.items(), key=lambda t: (len(t[0]), t[0]))
    names = {}
    out = []
    for m, c in ms:
        if not m:
            out.append('K' if abs(c) > 1000 else str(c))
            continue
        vs = []
        for v in m:
            if v not in names: names[v] = chr(ord('a')+len(names))
            vs.append(names[v])
        cc = 'K' if abs(c) > 1000 else str(c)
        out.append(cc + '*' + '*'.join(sorted(vs)))
    return ' '.join(out)

def main():
    sc = collections.Counter()
    for p in facs: sc[shape(p)] += 1
    print("distinct factor shapes:", len(sc))
    for s, n in sc.most_common(40):
        print("  %6d  %s" % (n, s))
    # sample of the rare ones
    print("--- rare shapes")
    rare = [s for s, n in sc.items() if n < 50]
    for s in rare[:40]:
        print("  %6d  %s" % (sc[s], s))

if __name__ == '__main__':
    main()
