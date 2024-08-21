$TITLE Riplito

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

model ripley /all/ ;

parameter rep ; 

solve ripley using lp maximizing z ; 
rep(i,"bau") = X.l(i) ; 

* enable the switch and re-solve
sw_combo = 1 ; 
solve ripley using lp maximizing z ; 
rep(i,"combo") = X.l(i) ; 

sw_combo = 0 ; 


c("hot_dogs") = 2 * c("hot_dogs") ; 
h("hot_dogs") = 2 * h("hot_dogs") ; 
solve ripley using lp maximizing z ; 
rep(i,"hungry") = X.l(i) ; 


execute_unload 'alldata.gdx' ; 