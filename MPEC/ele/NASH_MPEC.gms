$title A non-cooperative Game: Nash and Stackelberg versions (NASH,SEQ=269)

$onText
The original formulation is as a non-cooperative game;
a Nash equilibrium is sought.  The Nash game is then extended to a
Stackelberg or leader-follower game.  The Stackelberg game differs
from the original Nash game in that one firm (the leader)
now anticipates the prices set by all the other firms (the followers) and
sets its own prices based on this information.


References:

F.H. Murphy, H.D. Sherali, and A.L. Soyster, "A mathematical
programming approach for determining oligopolistic market
equilibrium", Mathematical Programming 24 (1982) 92-106

P.T. Harker, "Accelerating the Convergence of the Diagonalization and
Projection Algorithms for Finite-Dimensional Variational Inequalities",
Mathematical Programming 41 (1988) 29-59.

S.P. Dirkse and M.C. Ferris, "MCPLIB: A Collection of Nonlinear Mixed
Complementarity Problems", Optimization Methods in Software 5 (1995), 319-345.


Together with a production quantity vector q(I), the data below define:

0. A total quantity TQ = sum{I,q(I)} being produced

1. An inverse demand function p(TQ), i.e. the unit price at which consumers
   will demand (and actually purchase) a quantity TQ.  The inverse demand
   function p(TQ) is derived from the demand curve

     TQ(p) = dbar*p**(-gamma)

   yielding

     p(TQ) = dbar**(1/gamma)*TQ**(-1/gamma)

2. The total production cost TC_i(q_i) for firm i (not the per-unit cost)
   If we equate the marginal cost with price and assume each firm's supply
   curve is a constant elasticity function q = L*p**beta for some constant L
   and supply elasticity beta, we get a price function

     p = (q/L)**(1/beta)

   If we add a constant c to the price above (a fixed per/unit cost), the
   following function TC has the desired marginal cost (i.e. price):

     TC(I) = c(I)*q(I) +  beta(I)/(1+beta(I))*L(I)**(1/beta(I))
                                   *q(I)**((1+beta(I))/beta(I))


Individual firms that behave in a Nash manner
(i.e. they assume the other firms' decisions are fixed) optimize their
production q(I) to maximize profit:

  profit(i) = q(I)*p(TQ) - TC(I)

Keywords: mixed complementarity problem, equilibrium constraints, Nash game,
          Stackelberg game, oligopolistic market equilibrium, general equilibrium
          model
$offText

Set
   I         'firms' / 1*10 /
   leader(I) 'at most one leader firm';

Table IData(I,*) 'data for the firms'
        c   L  beta
    1   5  10  1.2
    2   3  10  1
    3   8  10   .9
    4   5  10   .6
    5   1  10  1.5
    6   3  10  1
    7   7  10   .7
    8   4  10  1.1
    9   6  10   .95
   10   3  10   .75;

Scalar
   gamma   'demand elasticity'               /  1.2 /
   dbar    'reference or unit-price demand'  / 5000 /;

Parameter
   c(I)    'unit production cost'
   L(I)    'supply curve scale factor'
   beta(I) 'supply curve elasticity'
   s(I);

c(I)    = IData(I,'c');
L(I)    = IData(I,'L');
beta(I) = IData(I,'beta');
s(I)    = (beta(I)/(1 + beta(I)))*L(I)**(1/beta(I));

Alias (I,J,K);

Positive Variable
   q(I) 'production vector';

Variable
   p 'price'
   z 'leader profit';

Equation
   dNegProfit(I) 'd(-profit)/dq: dTC - p - q(I)*dp'
   obj           'leader profit = q(I)*p - TC'
   defp          'total production';

defp..
   p =e= (dbar/sum(J, q(J)))**(1/gamma);

dNegProfit(I)$[not leader(I)] ..
   c(I) + (L(I)*q(I))**(1/beta(I)) - p + q(I)*(1/gamma)*p/sum(J, q(J)) =g= 0;

obj..
   z =e= sum{i$leader(i), q(I)*p - q(I)*c(I) - s(I)*q(I)**((1 + beta(I))/beta(I))};

Model
   nashmcp  / dNegProfit.q, defp.p      /
   nashmpec / dNegProfit.q, defp.p, obj /;

* allow some domain violations in DNLP subsolvers
option domLim = 50;
$ifi %gams.mcp% == nlpec option dnlp = conopt;

Parameter report;

p.l       =  1;
q.l(I)    =  1;
leader(I) = no;
solve nashmcp using mcp;

* compute player profits in the Nash case,
* i.e. when all players assume nothing about the other players' strategies
report(I,'nash') = q.l(I)*p.l - c(I)*q.l(I) - s(I)*q.l(I)**((1 + beta(I))/beta(I));
report('price','nash') = p.l;

leader(K) = no;
loop(K,
   leader(K) = yes;
   solve nashmpec using mpec max z;

*  compute player profits in the Stackelberg (leader-follower) case
   report(I,K) = q.l(I)*p.l - c(I)*q.l(I) - s(I)*q.l(I)**((1 + beta(I))/beta(I));
   report('price',K) = p.l;
   leader(K) = no;
);

Parameter reprel 'percentage change';
reprel(I,J)       = 100*(report(I,J) - report(I,'nash'))/report(I,'nash');
reprel('price',J) = 100*(report('price',J) - report('price','nash'))/report('price','nash');

display report, reprel;