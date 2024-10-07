
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

objfn.. W =e= beta * (
            ((Qd("pizza"))**alpha("pizza"))
            * ((Qd("beer"))**alpha("beer"))
            ); 

sup_dem(i).. Qs(i) =e= Qd(i) ; 

income_limit.. income =g= sum(i, pbar(i) * Qd(i)) ; 

model welfy /all/ ; 

Qd.lo(i) = 0.01 ; 
Qs.lo(i) = 0.01 ; 

solve welfy using nlp maximizing W ; 

execute_unload 'alldata.gdx' ; 
