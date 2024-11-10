$Title Simple Hourly Dispatch LP
* Maxwell Brown

$if not set sw_season $setglobal sw_season "day"

*Policy Switches
Scalar CAPTrade "Switch to turn on [1] or off [0] a carbon cap with trading" /0/;
Scalar TPSTrade "Switch to turn on [1] or off [0] a carbon TPS with trading" /0/;

*=======================
*  Begin sets
*=======================

set h "hour" /h1*h24/;
set k "season" /day, summer/;

*table d(h,k) "demand by season, MWH"
Table dem_in(h,k) "reference demand by hour and season (MWh)"
$include refdem.inc
;

parameter d(h) "demand by hour - used in the model (MWh)";
d(h) = dem_in(h,"%sw_season%");

*!!!!
set a 'scen' /a1*a4/ ; 
set dem_var /p,q/ ; 


table ref_pq(a,h,dem_var)
$ondelim
$include dem_scen.csv
$offdelim
; 
*!!!!

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

Table plantdata(s,f,pid,pc)
$include plantdata.inc
;

*set the generation feasibility set such that it is enabled
*when the indicated s/f/pid/h combination has capacity
genfeas(s,f,pid,h)$(plantdata(s,f,pid,"CAP")>0) = yes;
*can't use solar when it's dark out
genfeas(s,"sun",pid,h)$((ord(h)<8) or (ord(h)>18)) = NO;

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
*taken from EPA's model legislation
*see results browser documentation for table

*Policy Parameters
*Unless otherwise defined, the Psi parameters are relative to the base case
*note that we cannot be too stringent w/o elastic demand, safety valve credits, or some capacity expansion capabilities
Parameter  psi(s) "Percentage reduction of average CI or total carbon" 
                /CO 10, WY 7.5/;

psi(s) = psi(s) / 100;

*These start out as large values just in case the switch is activated
*without being properly calibrated to the reference case's reduction
*thus making it a non-binding constraint
Parameter ebar_cap(s) "Reference Aggregate Carbon level at which the standard is met";
Parameter ebar_rate(s) "Reference average carobn intensity at which the standard is met";
ebar_cap(s) = 1e10;
ebar_rate(s) = 1e10;

*=======================
*  End data
*=======================


*=======================
*  Begin Model
*=======================

equation 
eq_capcon(a,s,f,pid,h) "generation cannot exceed capacity",
eq_demcon(a,h) "generation must equal demand",
eq_demcon_mcp(a,h) 
zpc_dual;

positive variable X(a,s,f,pid,h) "generation in MWH",
lambda_k, lambda_d, lambda_tps, lambda_cap(s) ; 

variable eta_est ; 

variable Z "objective function value ($s)";

*!!!!
equation eq_obj_mpec ; 
eq_obj_mpec.. Z =e=  1e-4 * 
                sum((a,h),power(ref_pq(a,h,"p") - lambda_d(a,h),2) );
;


eq_capcon(a,s,f,pid,h)$genfeas(s,f,pid,h).. 
        capfac(f) * plantdata(s,f,pid,"cap") =g= X(a,s,f,pid,h);

scalar eta /-.5/ ; 
eq_demcon_mcp(a,h).. 
                sum((s,f,pid)$genfeas(s,f,pid,h),X(a,s,f,pid,h)) 
*!!!! note - a1 is our reference scenario
                =e= 
                ref_pq("a1",h,"q") * (ref_pq(a,h,"p") / ref_pq("a1",h,"p")) ** (eta); 

equation eq_demcon_mpec ; 

eq_demcon_mpec(a,h).. 
                sum((s,f,pid)$genfeas(s,f,pid,h),X(a,s,f,pid,h)) 
*!!!! note - a1 is our reference scenario
                =e= 
                ref_pq("a1",h,"q") * (ref_pq(a,h,"p") / ref_pq("a1",h,"p")) ** (eta_est); 


zpc_dual(a,s,f,pid,h)$genfeas(s,f,pid,h)..
             (plantdata(s,f,pid,"onm") + plantdata(s,f,pid,"hr") * fcost(f))
             +  lambda_k(a,s,f,pid,h)
             =g= 
             lambda_d(a,h) ; 


*!!!!
*lambda_k.l(a,s,f,pid,h)$genfeas(s,f,pid,h) = eq_capcon.m(s,f,pid,h) ; 
lambda_d.l(a,h) = ref_pq(a,h,"p"); 
lambda_d.lo(a,h) = 0.5 * ref_pq(a,h,"p"); 
lambda_d.up(a,h) = 2 * ref_pq(a,h,"p"); 




model mcp_ele 
/
eq_demcon_mcp.lambda_d,
eq_capcon.lambda_k,
zpc_dual.x
/ ; 

model mpec_ele 
/
eq_obj_mpec,
eq_demcon_mpec.lambda_d,
eq_capcon.lambda_k,
zpc_dual.x
/ ; 


solve mcp_ele using mcp ; 

eta_est.lo = -5 ; 
eta_est.up = -0.001 ; 


solve mpec_ele using mpec minimizing z; 

execute_unload 'refdata.gdx' ; 

$exit

ebar_cap(s) = sum((f,pid,h)$genfeas(s,f,pid,h), plantdata(s,f,pid,"hr") * emit(f) * X.l(s,f,pid,h) );

ebar_rate(s) = sum((f,pid,h)$genfeas(s,f,pid,h), plantdata(s,f,pid,"hr") * emit(f) * X.l(s,f,pid,h) )
                / sum((f,pid,h)$genfeas(s,f,pid,h), X.l(s,f,pid,h) )
;

lambda_d.lo(h) = 1e-3 ; 

mcp_ele.iterlim = 1000000 ; 

sw_elas = 1 ; 
swtps(s) = yes ; 
swcap(s) = no ;

solve mcp_ele using mcp ; 





