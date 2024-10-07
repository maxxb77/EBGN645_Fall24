
set i /pizza, beer/ ; 

parameter c0(i), d0(i); 

d0(i) = 2 ;

parameter pbar /pizza 2, beer 1/ ; 
parameter qbar(i) /pizza 3, beer 1/ ; 

c0(i) = qbar(i) - d0(i) * pbar(i) ; 

scalar income ; 

income = sum(i,qbar(i) * pbar(i)) ; 

parameter c(i), d(i) ; 

c(i) = -(c0(i)/d0(i)) ;
d(i) = 1/d0(i) ; 


positive variables 
Qd(i) "quantity demanded"
Qs(i) "quantity supplied"
;

variable W "social welfare" ; 

equation objfn, sup_dem(i), income_limit ; 

scalar beta /1/ ; 
parameter alpha(i); 
alias(i,ii) ; 

alpha(i) = (pbar(i) * qbar(i))/(sum(ii,pbar(ii) * qbar(ii))) ; 
beta = sum(ii,pbar(ii) * qbar(ii)) ; 
scalar rho /0.5/ ; 

objfn.. W =e= beta * 
            (sum(i,alpha(i) * (qd(i)/qbar(i)) **rho)
            )**(1/rho); 

sup_dem(i).. Qs(i) =e= Qd(i) ; 

parameter priceadder(i) ; 
priceadder(i) = 0 ; 

income_limit.. income =g= sum(i, (1+priceadder(i)) * pbar(i) * Qd(i)) ; 

model welfy /all/ ; 

Qd.lo(i) = 0.01 ; 
Qs.lo(i) = 0.01 ; 

parameter report ; 

solve welfy using nlp maximizing W ; 
report("ref","consumption",i) = qd.l(i) ; 
report("ref","welfare","welfare") = w.l ; 

priceadder("pizza") = 0.1 ; 
solve welfy using nlp maximizing W ; 
report("tax","consumption",i) = qd.l(i) ; 
report("tax","welfare","welfare") = w.l ; 


execute_unload 'alldata.gdx' ; 
