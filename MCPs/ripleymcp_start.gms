$TITLE Riplitoni

$ontext
this is a really fun class and i love hamburgers
$offtext

* this is my set
set i "items produced" 
    / hot_dogs, hamburgers, frenchies / ; 

parameter 
r(i) "--$s/unit-- revenue per unit of production"
/
hamburgers 1.5,
hot_dogs 0.75,
frenchies 0.25
/,
c(i) "--$s/unit-- cost per unit of production"
/
hamburgers 0.75,
hot_dogs 0.2,
frenchies 0.1
/
h(i) "--hours/unit-- hours required per unit"
/
hamburgers 0.5,
hot_dogs 0.2,
frenchies 0.1
/ ; 

scalar hbar "total hours worked in a week" /40/ ; 

positive variable x(i) "--units-- production of different products" ; 
variable z "--$s-- total profits" ; 

equation objfn "calculation of our maximization target", 
eq_hourlimit "can only work 40 hours in a week",
eq_combo "french fries production needs to be equal to or greater than hot dog production" ; 

scalar sw_combo "boolean switch to turn on and off the combo constraint" /0/ ; 

objfn.. z =e= sum(i, (r(i) - c(i)) * X(i) ) ;

eq_hourlimit.. hbar =g= sum(i,h(i) * X(i) ) ; 

eq_combo$sw_combo.. X("frenchies") =g= X("hot_dogs") ; 

model rip_lp /all/ ;

* for converting LP/NLP into MCP...
*1 solve the linear program
sw_combo = 1 ; 
solve rip_lp maximizing Z using lp ;


*2 set level values to marginals off constraints
* noting levels may be off in sign given problem direction
positive variable lambda_hour, lambda_combo ;
lambda_hour.l = abs(eq_hourlimit.m) ; 
lambda_combo.l = abs(eq_combo.m) ; 


* set up zero profit condition
equation zpc_x(i) ;

zpc_x(i)..
 lambda_hour * h(i) 
 + lambda_combo$[sameas(i,"hot_dogs")$sw_combo]
 =g=
 (r(i)-c(i))
 + lambda_combo$[sameas(i,"frenchies")$sw_combo]
;

*!!! be careful on switch settings and dual value inclusion

* set up complementarity relationships
model ripley_mcp
/
eq_hourlimit.lambda_hour,
eq_combo.lambda_combo,
zpc_x.X
/;

* set iteration limit to zero and solve
* making sure we have a zero residual
ripley_mcp.iterlim = 0 ; 

solve ripley_mcp using mcp ; 



