
set i /pizza, beer/ ; 

parameter c0(i), d0(i); 

d0(i) = 2 ;

scalar pbar /1/ ; 
parameter qbar(i) /pizza 3, beer 1/ ; 

c0(i) = qbar(i) - d0(i) * pbar ; 

parameter c(i), d(i) ; 

c(i) = -(c0(i)/d0(i)) ;
d(i) = 1/d0(i) ; 


positive variables 
Qd(i) "quantity demanded"
Qs(i) "quantity supplied"
;

variable W "social welfare" ; 

equation objfn, sup_dem(i) ; 

scalar income ; 

income = sum(i,pbar * qbar(i)) ; 

scalar beta /4/ ; 
parameter theta(i) /pizza 0.75, beer 0.25/ ; 

parameter priceadj(i);
priceadj(i) = 0 ; 


scalar rho /2/ ; 

objfn.. W =e= beta * (
                sum(i, theta(i) * (qd(i)/qbar(i))**rho )**(1/rho)
*            ((Qd("pizza")/qbar("pizza"))**alpha("pizza"))
*            * ((Qd("beer")/qbar("beer"))**alpha("beer"))
            )
*            - sum(i, (c(i) * Qs(i) + d(i) * Qs(i) * Qs(i) / 2))
; 

equation 
income_constraint ; 

income_constraint.. income =g= sum(i,(1+ priceadj(i)) * pbar * qd(i)) ; 

sup_dem(i).. Qd(i) =e= Qs(i) ; 

model welfy /all/ ; 

Qd.lo(i) = 0.01 ; 
Qs.lo(i) = 0.01 ; 

parameter report;

solve welfy using nlp maximizing W ; 
report('ref',i,"consumption") = qd.l(i) ; 
report('ref',"na","utility") = w.l ; 

priceadj("pizza") = 0.1 ; 

solve welfy using nlp maximizing W ; 
report('tax',i,"consumption") = qd.l(i) ; 
report('tax',"na","utility") = w.l ; 

income = 1.5 * income ; 
priceadj(i) = 0 ; 

solve welfy using nlp maximizing W ; 
report('income',i,"consumption") = qd.l(i) ; 
report('income',"na","utility") = w.l ; 


execute_unload 'alldata.gdx' ; 
