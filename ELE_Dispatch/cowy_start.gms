$Title Simple Hourly Dispatch LP
* Maxwell Brown

option seed=7

$if not set sw_season $setglobal sw_season "day"
$if not set sw_co_reduction $setglobal sw_co_reduction 10
$if not set sw_wy_reduction $setglobal sw_wy_reduction 7.5

* enable elastic demand in the counterfactuals?
$if not set sw_elas $setglobal sw_elas 0



*=======================
*  Begin sets
*=======================

set h "hour" /h1*h24/;
set dayhours(h) /h8*h20/; 

set k "season" /day, summer/;

set s "states"
/
  CO "Colorado",
  WY "Wyoming"
/;

*State-specific policy switches
Set SwTPS(s)  "Switch to enable or disable a TPS for state s";
Set SwCAP(s)  "Switch to enable or disable a carbon cap for state s";

*disable these for now
swtps(s) = no;
swcap(s) = no;

set f "fuels, generation technology"
        /
        bit "bituminous coal"
        dfo "distillate fuel oil"
        ng  "natural gas"
        sub "subbituminous coal"
        sun "solar"
        wat "water"
        wnd "wind"
        /;

set pc  "plant characteristics" /cap, hr, onm/;
set pid "plant id" /1*328/;
set genfeas(s,f,pid,h) "generation feasibility set, determines when a plant can generate";

*=======================
*  End sets
*=======================

*=======================
*  Begin data
*=======================

* some rudimentary data...
parameter capfac(f) "capacity availability by fuel type"
         /
         bit 0.8,
         dfo 0.65,
         ng  0.85,
         sub 0.8,
         sun 0.5,
         wat 0.55,
         wnd 0.45
         /;
* offhand guesses 

parameter cf(f,h) 'capacity factor by technology and hour' ; 

* start by broadcasting values
cf(f,h) = capfac(f) ; 

* remove those not possible
cf("sun",h)$(not dayhours(h)) = 0 ; 


parameter fcost(f) "costs by fuel ($ / MMBTU)"
         /
         bit  2.9,
         dfo  11.1,
         ng   2.97,
         sub  2.9
         /;
*estimates from EIA, taken from short term energy outlook

parameter emit(f) "lbs co2 per mmbtu for fuel burning"
         /
         bit  202,
         dfo  161,
         ng   117,
         sub  209
         /;

parameter 
  psi_state(s), 
  psi_trade, 
  phi_state(s), 
  phi_trade ; 

* set all to zero to start
psi_state(s) = 0 ;  
psi_trade = 0 ; 
phi_state(s) = 0 ; 
phi_trade = 0 ; 


table d_in(h,k) 
$include refdem.inc
; 

parameter d(h) "demand by hour (MW)" ; 

d(h) = d_in(h,"%Sw_Season%") ; 

table plantdata(s,f,pid,pc)
$include plantdata.inc
;

set v(s,f,pid,h) "allowable combinations for technologies over all indices"; 
*populate combinations..
* - if you have a capacity factor (specific to solar)
* - if you have capacity
v(s,f,pid,h)$[cf(f,h)$plantdata(s,f,pid,"cap")] = yes ; 

scalar sw_elas /0/ ; 

set b 'demand bin' /1*100/ ; 
alias(b,bb) ; 
scalar ele_elas /-0.5/ ; 

parameter pbar(b,h) "value of demand for each bin and hour ($ / mwh)",
          qbar(b,h) "quantity available for each bin and hour (mwh)"
;

pbar(b,h) = 0 ; 
qbar(b,h) = d(h) / 50 ; 


execute_unload 'demand_check.gdx' ; 

positive variables 
X(s,f,pid,h) "generation by unit",
Q(b,h) ; 

variable 
Z "total cost - target of our objective" ; 


equation 
eq_objfn, eq_caplimit(s,f,pid,h), eq_supply_demand(h), eq_qbar_limit(b,h) ; 

eq_objfn.. 
  Z =e= sum((s,f,pid,h)$v(s,f,pid,h), X(s,f,pid,h) * 
    (plantdata(s,f,pid,"onm") + plantdata(s,f,pid,"hr") * fcost(f)) )
    - sum((b,h),pbar(b,h) * Q(b,h))$sw_elas
; 

eq_qbar_limit(b,h)$sw_elas.. qbar(b,h) =g= Q(b,h) ; 

eq_caplimit(s,f,pid,h)$v(s,f,pid,h)..
  cf(f,h) * plantdata(s,f,pid,"cap") =g= X(s,f,pid,h) 
; 

eq_supply_demand(h)..
  sum((s,f,pid)$v(s,f,pid,h), X(s,f,pid,h)) =e= 
  d(h)$(not sw_elas)
  + sum(b,Q(b,h))$sw_elas ; 

*Policy Switches
scalar
  cap_state "impose cap on state level" /0/
  cap_trade "impose cap with trading"   /0/
  tps_state "impose tps on state level" /0/
  tps_trade "impose tps with trading"   /0/
;

* policy equations
equation 
eq_cap_state(s)
eq_cap_trade
eq_tps_state(s)
eq_tps_trade
;

parameter phi, phi_state, psi, psi_state ; 

phi = 0 ; 
phi_state(s) = 0 ; 
psi = 0 ; 
psi_state(s)  = 0 ; 

eq_cap_state(s)$cap_state.. 
  psi_state(s) 
  
  =g=

  sum((f,pid,h)$v(s,f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X(s,f,pid,h)) 
;


eq_cap_trade$cap_trade.. 

  psi
  
  =g=

  sum((s,f,pid,h)$v(s,f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X(s,f,pid,h)) 
;


eq_tps_state(s)$tps_state.. 

  phi_state(s) * sum((f,pid,h)$v(s,f,pid,h), X(s,f,pid,h)) 
  
  =g=

  sum((f,pid,h)$v(s,f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X(s,f,pid,h)) 


;

eq_tps_trade$tps_trade.. 

  phi * sum((s,f,pid,h)$v(s,f,pid,h), X(s,f,pid,h)) 
  
  =g=

  sum((s,f,pid,h)$v(s,f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X(s,f,pid,h)) 


;

model ele /all/ ; 

parameter rep_gen ; 
parameter rep_gen_day, rep_cost ; 

solve ele using lp minimizing z ; 

* compute the values along the demand curve
pbar(b,h) = eq_supply_demand.m(h) * 
  (sum(bb$(bb.val<=b.val),qbar(b,h)) / d(h)) ** ele_elas;

* enable elasticity
sw_elas = %sw_elas% ; 

* solve our reference case
solve ele using lp minimizing z ; 
rep_gen(s,f,pid,h,"bau") = X.l(s,f,pid,h) ; 
rep_gen_day(s,f,"bau") = sum((pid,h),X.l(s,f,pid,h)) ; 
rep_cost("bau") = z.l ;

scalar 
co_reduction /%sw_co_reduction%/ 
wy_reduction /%sw_wy_reduction%/ ; 

parameter state_red(s) ; 
state_red("co") = 1-co_reduction/100 ; 
state_red("wy") = 1-wy_reduction/100 ; 

psi_state(s) = 
    state_red(s) * sum((f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X.l(s,f,pid,h) )
    ;

cap_state = 1 ; 

solve ele using lp minimizing z ; 
rep_gen(s,f,pid,h,"state_cap") = X.l(s,f,pid,h) ; 
rep_gen_day(s,f,"state_cap") = sum((pid,h),X.l(s,f,pid,h)) ; 
rep_cost("state_cap") = z.l ;


* grab state-level emissions rates 
phi_state(s) =
  sum((f,pid,h)$v(s,f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X.l(s,f,pid,h) )
  / sum((f,pid,h)$v(s,f,pid,h), X.l(s,f,pid,h) )
  ;


* disable state cap
* enable trade cap
cap_state = 0 ; 
cap_trade = 1 ; 
psi = sum(s,psi_state(s)) ; 

solve ele using lp minimizing z ; 
rep_gen(s,f,pid,h,"trade_cap") = X.l(s,f,pid,h) ; 
rep_gen_day(s,f,"trade_cap") = sum((pid,h),X.l(s,f,pid,h)) ; 
rep_cost("trade_cap") = z.l ;

phi =
  sum((s,f,pid,h)$v(s,f,pid,h), emit(f) * plantdata(s,f,pid,"hr") * X.l(s,f,pid,h) )
  / sum((s,f,pid,h)$v(s,f,pid,h), X.l(s,f,pid,h) )
  ;

*disable tradable cap
cap_trade = 0 ; 

*enable state-level emissions standards
tps_state = 1 ; 

solve ele using lp minimizing z ; 
rep_gen(s,f,pid,h,"tps_state") = X.l(s,f,pid,h) ; 
rep_gen_day(s,f,"tps_state") = sum((pid,h),X.l(s,f,pid,h)) ; 
rep_cost("tps_state") = z.l ;

*disable state-level tps
tps_state = 0 ;

* enable tradable performance standard for all states 
tps_trade = 1 ;

solve ele using lp minimizing z ; 
rep_gen(s,f,pid,h,"tps_trade") = X.l(s,f,pid,h) ; 
rep_gen_day(s,f,"tps_trade") = sum((pid,h),X.l(s,f,pid,h)) ; 
rep_cost("tps_trade") = z.l ;


execute_unload 'cowy_%elas%.gdx' ; 