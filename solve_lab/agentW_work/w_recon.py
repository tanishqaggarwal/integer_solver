"""W: reconcile U's five discount equations {2554,6816,8124,8680,9421} with MY essential-row
family {2554,6816,8124,9123,9421,S} and my round-1 32-way trade prices.  Is 'S' eq8680?"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
d = model.get(); EQ = d['eq_terms']; A = d['atom_src']
SATOM = 37887
eqs_with_S = [i for i, (m, sq, tl) in enumerate(EQ) if any(a == SATOM for c, a in tl)]
print('atom a%d (my "S" row) appears in equations: %s' % (SATOM, eqs_with_S))
print('  is 8680 among them? ', 8680 in eqs_with_S)
m, sq, tl = EQ[8680]
print('  eq8680 has %d atoms, square=%s ; coefficient on a37887: %s'
      % (len(tl), sq, [c for c, a in tl if a == SATOM]))
MINE_ESS = ['2554', '6816', '8124', '9123', '9421', 'S']
MINE_PRICES = [2554, 6816, 8124, 9123, 9421, 8680]     # round-1 audit, the 32-way trade
U_FIVE = [2554, 6816, 8124, 8680, 9421]
print()
print('my essential rows        :', MINE_ESS)
print('my round-1 trade prices  :', sorted(MINE_PRICES))
print("U's five discount eqs    :", sorted(U_FIVE))
print("U's five  subset of  my prices? ", set(U_FIVE) <= set(MINE_PRICES),
      ' ; my prices minus U =', sorted(set(MINE_PRICES) - set(U_FIVE)))
print()
# is 8680 itself one of the 205 rows of the K+ system, and is it SAT or FAIL?
import w_setup2 as S
print('rows in the K+ system: %d ; is 8680 a row? %s ; 9123? %s'
      % (len(S.names), 8680 in S.names, 9123 in S.names))
print('  8680 in SAT? %s   in FAIL? %s' % (8680 in S.SAT, 8680 in S.FAIL))
print('  dropped-as-nonlinear rows:', S.dropped)
print('  is 8680 among the dropped?', 8680 in S.dropped)
