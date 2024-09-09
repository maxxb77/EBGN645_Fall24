$TITLE Ursula

set i / pizza, beer / ;

parameters 
u(i) "utility per item consumed (utils/item)"
/
pizza 2
beer 5
/, 

cost(i) "cost per item ($s/item)"
/
pizza 2
beer 4
/ ;

scalar b "budget ($s)" /20/ ;

positive variable X(i) "consumption of each item (items)" ; 
variable utility "total utility (utils)" ; 

equations
eq_objfn "calculation of total utility", 
eq_budget "budget limit" ; 

eq_objfn.. utility =e= sum( i , u(i) * X(i)) ; 

eq_budget.. b =g= sum(i, cost(i) * X(i)) ; 

model ursula /all/ ; 

solve ursula maximizing utility using lp ; 

execute_unload 'ursula.gdx' ; 

