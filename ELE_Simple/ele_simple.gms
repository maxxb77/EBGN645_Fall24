
set i /coal, gas, wind/,  
    r /r1, r2/,
    h /1*24/ ; 

alias(r,rr) ; 
alias(h,hh) ; 

* have normally had simple examples like this...
parameter c(i)
/
Coal	20
Gas	15
Wind	1
/ ;

set vre(i) /wind/ ; 

* can also store data in a .gms file:
parameter cf_wind(h) "wind capacity factor"
* include an already-formatted-for-gams table
$include wind.gms
;

parameter cf(i,h) ; 

cf(i,h)$[vre(i)] = cf_wind(h) ; 
cf(i,h)$[not vre(i)] = 1 ; 

* load in a csv
table demand(h,r)
* tell gams to see commas as separators of data
$ondelim
$include demand.csv
$offdelim
;

table capacity(i,r)
$ondelim
$include capacity.csv
$offdelim
;

set v(i,r,h) "valid combination of i/r/h" ; 
v(i,r,h)$[capacity(i,r)$cf(i,h)] = yes ; 

scalar transmission_capacity "capacity of transmisison between both regions" /5/ ; 
scalar transmission_loss "loss rate from one region to another" /0.05/ ; 

positive variable
G(i,r,h) "generation (MWh)"
S(r,rr) "transmission from r to rr" 
*DEMAND_SLACK(r,h) "slack variable for demand"
; 

variable 
Z "cost - target of our optimization"
;

Equations
eq_objfn "objective function"
eq_caplimit(i,r,h) "generation cannot exceed capacity"
eq_transmission_limit "transmission cannot exceed transmission capacity"
eq_demand(r,h) "supply (generation [plus] net transmission) must equal demand"
;

eq_objfn..
    Z =e= sum((i,r,h)$v(i,r,h), c(i) * G(i,r,h) ) 
*    + sum((r,h), DEMAND_SLACK(r,h) * 1e5)
; 

eq_caplimit(i,r,h)$v(i,r,h)..
    cf(i,h) * capacity(i,r) =g= G(i,r,h) ;

eq_transmission_limit..
    transmission_capacity =g= sum((r,rr)$[not sameas(r,rr)], S(r,rr) + S(rr,r) ) ;

eq_demand(r,h)..
* domestic generation
*    DEMAND_SLACK(r,h)

    + sum((i)$v(i,r,h), G(i,r,h) ) 
* net transmission
    + sum(rr$[not sameas(r,rr)], (1-transmission_loss) * S(rr,r) - S(r,rr) )

    =g= 

    demand(h,r)
;

parameter delta(i)
/
coal 0.02, 
gas 0.1
/;

equation eq_unitcommitment_lower(i,r,h,hh), eq_unitcommitment_upper(i,r,h,hh)  ; 

set eqlimit(i,r,h,hh) ; 
eqlimit(i,r,h,hh)$[v(i,r,h)$v(i,r,hh)$(not sameas(h,hh))$delta(i)] = yes ; 

eq_unitcommitment_lower(i,r,h,hh)$eqlimit(i,r,h,hh)..
    G(i,r,h)
    =g=
    (1-delta(i)) * G(i,r,hh)
;

eq_unitcommitment_upper(i,r,h,hh)$eqlimit(i,r,h,hh)..


    (1+delta(i)) * G(i,r,hh)
    =g=
    G(i,r,h)
;


model elmo /all/ ; 

solve elmo using lp minimizing z ; 

* is wind generation being curtailed?
* basically: is the output less than the capacity?

parameter rep_wind ; 
rep_wind(r,h) = g.l("wind",r,h) / (cf("wind",h) * capacity("wind",r)) ; 

execute_unload 'data.gdx' ; 