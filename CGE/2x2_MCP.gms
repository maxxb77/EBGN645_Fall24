* A 2 sector, 2 factor CGE model based on the SAM presented in lecture.

$ontext
                  Production Sectors          Consumers
   Markets   |    X       Y        W    |       CONS
   ------------------------------------------------------
        PX   |  100             -100    |
        PY   |          100     -100    |
        PW   |                   200    |       -200
        PL   |  -40     -60             |        100
        PK   |  -60     -40             |        100
   ------------------------------------------------------
$offtext

SCALARS
	tx	tax rate in the X sector	/0/
	lendow	labor endowment multiplier	/1/
;

POSITIVE VARIABLES
	 PX	price of the X good
	 PY	price of the Y good
	 PW	price index for consumer
	 PL	price of labor
	 PK	price of capital

	 X	activity level in X sector
	 Y	activity level in Y sector
	 W	utility index
	 CONS	consumer income level
;

EQUATIONS
	PRF_X	zero profit condition for X sector
	PRF_Y	zero profit condition for Y sector
	PRF_W	zero profit condition for utility

	MKT_X	market clearance for X
	MKT_Y	market clearance for Y
	MKT_L	market clearance for L
	MKT_K	market clearance for K
	MKT_W	market clearance for utility good

	I_CONS	income balance condition for consumer
;

PRF_X..		PL**0.4 * PK**0.6 * (1 + tx) =g= PX;

PRF_Y..		PL**0.6 * PK**0.4 =g= PY;

PRF_W..		PX**0.5 * PY**0.5 =g= PW;

MKT_X..		100*X =g= 100 * W * PX**0.5 * PY**0.5 / PX;

MKT_Y..		100*Y =g= 100 * W * PX**0.5 * PY**0.5 / PY;

MKT_W..		200*W =g= CONS / PW;

MKT_L..		100*lendow =g= 40 * X * PL**0.4 * PK**0.6 / PL
			       + 60 * Y * PL**0.6 * PK**0.4 / PL ;

MKT_K..		100 =g= 60 * X * PL**0.4 * PK**0.6 / PK
			+ 40 * Y * PL**0.6 * PK**0.4 / PK ;

* Tax revenue is returned lump-sum to the household.
I_CONS..	CONS =g= 100*lendow*PL + 100*PK + tx * 100 * X * PL**0.4 * PK**0.6;

MODEL twobytwo /PRF_X.X, PRF_Y.Y, PRF_W.W, MKT_X.PX, MKT_Y.PY,
      	       MKT_W.PW, MKT_L.PL, MKT_K.PK, I_CONS.CONS/;

* Replicate the benchmark equilibrium.

** Load benchmark prices and quantities.
X.l = 1; Y.l = 1; W.l = 1; PX.l = 1; PY.l = 1; PK.l = 1; PW.l = 1;

** Make the wage rate the numeraire price.
PL.fx = 1;

CONS.l = 200;

twobytwo.iterlim = 0;
solve twobytwo using MCP;




* Counterfactuals:

*introduce a 50 pct tax in the X sector.

twobytwo.iterlim = 1000;

tx = 0.5;
solve twobytwo using MCP;


** Second, double the size of the labor endowment.

tx = 0;
lendow = 2;

solve twobytwo using MCP;
