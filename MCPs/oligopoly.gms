
Set i 'firms' / i1*i5 /;

* f (q) = c q + beta/(beta+1) l^(-1/beta) q^(beta+1)/beta
* f'(q) = c + (q/l)^(1/beta)

Table data(i,*) 'cost function data'
         c  qbar  beta
   i1   10  1     1.2
   i2    8  1     1.1
   i3    6  1     1.0
   i4    4  1     0.9
   i5    2  1     0.8;

Parameter
   c(i)
   qbar(i)
   beta(i);

c(i)    = data(i,"c");
qbar(i)  = data(i,"qbar");
beta(i) = data(i,"beta");

Positive Variable
   p    'price'
   q(i) 'supply';

Equation
   demand    'supply - demand balance'
   profit(i) 'Nash first order condition';

*!!! fill in MCP model here

scalar alpha /5/ ;

scalar elas /1.5/ ; 

demand.. sum(i,q(i)) =g= alpha * P**(-elas) ; 

profit(i)..
   c(i) + (q(i)/qbar(i)) ** (1/beta(i))
   =g=
   P - alpha * q(i) * P**elas
; 


model oligopoly 
/
demand.P,
profit.q
/;


* initial guess:
q.l(i) =  1;
p.l    = 1 ;

solve oligopoly using mcp;

execute_unload 'alldata_cournot.gdx' ; 