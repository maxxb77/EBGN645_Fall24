* 2x2 CGE model with taxes in the X sector.

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
	tx	ad-valorem tax rate on X good	/0/
	lendow	labor endowment multiplier	/1/
;


* Now begin the description of the MPSGE model.

$ontext

$MODEL:2x2

$SECTORS:
	X	! Activity level for the sector X
	Y	! Activity level for the sector Y
	W	! Activity level for utility good

$COMMODITIES:
	PX	! Price index for the market for X
	PY	! Price index for the market for Y
	PL	! Wage rate
	PK	! Rental price of capital
	PW	! cost of living index

$CONSUMERS:
	CONS	! income level for the household CONS

$PROD:X	s:1
	o:PX	q:100	p:1	
	i:PL	q:40	a:CONS	t:tx
	i:PK	q:60	a:CONS	t:tx

$PROD:Y	s:1
	o:PY	q:100	p:1
	i:PL	q:60
	i:PK	q:40

$PROD:W	s:1
	o:PW	q:200
	i:PX	q:100
	i:PY	q:100

$DEMAND:CONS
	d:PW	q:200
	e:PL	q:(100*lendow)
	e:PK	q:100
	e:PX	q:100	r:XQADJ


$offtext

$sysinclude mpsgeset 2x2


* Choose numeraire price.

PL.fx = 1;

$include 2x2.gen
solve 2x2 using mcp;


* Solve the tax counterfactual:

2x2.iterlim = 1000;

tx = 0.5;

$include 2x2.gen
solve 2x2 using mcp;


