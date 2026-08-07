# U0: the three-line arithmetic, recomputed from scratch.
N = 115792089237316195423570985008687907852837564279074904382605163141518161494337
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
print("N bits          :", N.bit_length())
print("N popcount      :", bin(N).count('1'), " zeros:", 256-bin(N).count('1'))
print("2^256-1 >= 2N ? :", (2**256-1) >= 2*N)
print("2N - (2^256-1)  :", 2*N-(2**256-1))
print("2^256 - N       :", 2**256-N)
print("2^129           :", 2**129)
print("2^256-N < 2^129 :", (2**256-N) < 2**129)
# so maskval(J) >= N  <=>  J contains all of {129..255}?  Prove both directions numerically.
full = (1<<256)-1
# necessity: omit any e>=129  => maskval <= full - 2^e <= full - 2^129 < N ?
print("full - 2^129 < N:", (full - (1<<129)) < N)
# sufficiency: J = {129..255} alone
m = sum(1<<e for e in range(129,256))
print("maskval({129..255}) >= N :", m >= N, m-N)
# the unconstrained representability witness (coordinator / K): j with bit_j=1,bit_{j+1}=0
js=[j for j in range(255) if (N>>j)&1 and not (N>>(j+1))&1]
print("candidate j's (first 5):", js[:5], "count", len(js))
j=js[0]
S=set(e for e in range(256) if (N>>e)&1)
A=set(S); A.discard(j); A.add(j+1)   # 2^j -> 2^{j+1} - 2^j
B={j}
assert not (A & B), "A,B not disjoint"
val=sum(1<<e for e in A)-sum(1<<e for e in B)
print("witness j=%d  |A|=%d |B|=%d  sum_A - sum_B == N :"%(j,len(A),len(B)), val==N)
# also check: does that A contain all of {129..255}?
print("A superset of {129..255} :", set(range(129,256)) <= A)
