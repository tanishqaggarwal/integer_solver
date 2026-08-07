import pickle, os, collections
from model import Model
HERE = os.path.dirname(os.path.abspath(__file__))
M = Model()
val = pickle.load(open(os.path.join(HERE, 'prop0.pkl'), 'rb'))
known = [v for v in range(len(val)) if val[v] is not None]
big = [(v, val[v]) for v in known if abs(val[v]) > 10**9]
print("known:", len(known), " big-valued:", len(big))
c = collections.Counter(val[v].bit_length() for v in known if val[v] is not None)
print("bitlength histogram:", sorted(c.items())[:20], "...", sorted(c.items())[-10:])
seen = collections.Counter(val[v] for v, _ in big)
print("\ndistinct big constants:", len(seen))
for k, n in seen.most_common(40):
    print(f"  x{n:4d}  bits={k.bit_length():4d}  {k}")
