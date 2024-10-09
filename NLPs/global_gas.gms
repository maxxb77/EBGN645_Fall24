set r "region" /us, row/,
    i "energy goods" /gas, oth/ 
;

parameter
    b0(r,i)     "supply curve slope"
    /
    us.gas 0.25
    row.gas 0.75
    us.oth 0.75
    row.oth 0.5
    /
    pbar(r,i)    "reference global price"
    / us.gas 3, us.oth 3.2, row.gas 3, row.oth 3.2/
    qbar_consumption(r,i)  "reference quantity by region"
    /
* us data from:
* https://www.americangeosciences.org/critical-issues/faq/what-are-major-sources-and-users-energy-united-states
* world data from:
* https://www.iea.org/reports/world-energy-balances-overview/world
    us.gas 28, us.oth 70
    row.gas 73, row.oth 215
    /
    qbar_production(r,i)
    /
    us.gas 32, us.oth 70
* gas production computed as row consumption [minus] us net exports
    row.gas 69, row.oth 215
    /
    rho(r)     "elasticity of substitution across energy products"
    /us 0.5, row 0.5 /
* !!!
* calibrated values
    a0(r,i)     "inverse supply curve intercept"
    budget(r)  "total budget for energy goods"
    theta(r,i) "reference share of income spent on good i by region r"
;
b0(r,i) = 0.9 ; 
a0(r,i) = qbar_production(r,i) - pbar(r,i) * b0(r,i) ;
budget(r) = sum(i,qbar_consumption(r,i) * pbar(r,i)) ; 
theta(r,i) = qbar_consumption(r,i) * pbar(r,i) / budget(r) ; 


parameter a(r,i), b(r,i) ; 

a(r,i) = -(a0(r,i)/b0(r,i)) ;
b(r,i) = 1/b0(r,i) ; 

variable Z 'objective function target' ; 

positive variable production(r,i), consumption(r,i) ; 

equation objfn; 

parameter gamma(r) ; 
alias(r,rr) ; 
gamma(r) = budget(r) / sum(rr,budget(rr))
parameter beta ;
beta = 0.9; 

objfn.. Z =e= sum(r,budget(r)) * sum(r, gamma(r) * ((sum(i, theta(r,i) * (consumption(r,i) / qbar_consumption(r,i))**rho(r) ) 
                                )**(1/rho(r)) )**beta)**(1/beta)
              - sum((r,i), a(r,i) * production(r,i) + (b(r,i) * production(r,i)**2) / 2 )
;

equation supply_demand(i) ; 

supply_demand(i)..
    sum(r,production(r,i)) =e= sum(r,consumption(r,i)) ; 

model gas /all/ ;

consumption.lo(r,i) = 0.01 ; 
production.lo(r,i) = 0.01 ; 

parameter report ; 

solve gas using nlp maximizing Z ; 

report("bau","consumption",r,i) = consumption.l(r,i) ; 
report("bau","production",r,i) = production.l(r,i) ; 
report("bau","welfare",r,"N/A") = budget(r) * (sum(i,theta(r,i) * (consumption.l(r,i)/qbar_consumption(r,i))**rho(r) ) )**(1/rho(r)) ;


* gas shock 1 - us supply
*budget(r) = 1.001 * budget(r) ; 
*production.up("us","gas") = production.l("us","gas") * 0.9 ;
*consumption.lo("us","oth") = production.l("us","oth") * 1.1 ;  
*budget("us") = 1.1 * budget("us") ;  

a("us","gas") = 0.7 * a("us","gas") ;


solve gas using nlp maximizing Z ; 

report("gs1","consumption",r,i) = consumption.l(r,i) ; 
report("gs1","production",r,i) = production.l(r,i) ; 
report("gs1","welfare",r,"N/A") = budget(r) * (sum(i,theta(r,i) * (consumption.l(r,i)/qbar_consumption(r,i))**rho(r) ) )**(1/rho(r)) ;



execute_unload 'gas.gdx' ; 

