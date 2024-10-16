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



* --- Binning --- *
* create a set of bins

set b /1*100/ ; 
alias(b,bb) ; 

scalar sw_elas /0/ ; 
scalar elasticity /-2/ ; 


parameter pbar(h), qbar(h) ;

pbar(h) = 5 ; 
qbar(h) = d(h) ; 

parameter mu(h,b), qlim_bin(h,b) ; 
qlim_bin(h,b) = d(h) / 50 ; 
mu(h,b) = pbar(h) * (sum(bb$[bb.val<=b.val], qlim_bin(h,bb) ) / qbar(h) ) ** elasticity; 

*=======================
*  End data
*=======================


*=======================
*  Begin Model
*=======================

equation 
eq_cost "objective function",
eq_capcon(s,f,pid,h) "generation cannot exceed capacity",
eq_demcon(h) "generation must equal demand",
eq_dembin_lim(b,h)
eq_carboncap(s) "(optional) total emissions cannot exceed cap",
eq_ratestandard(s) "(optional) emissions rate cannot exceed emissions standard";

positive variable X(s,f,pid,h) "generation in MWH", 
                  dem_b(h,b) "demand by bin";

variable Z "objective function value ($s)";

eq_cost.. Z =e=  
        sum((s,f,pid,h)$genfeas(s,f,pid,h), 
                X(s,f,pid,h) * (plantdata(s,f,pid,"onm")
                                + plantdata(s,f,pid,"hr") * fcost(f)
                               )
           )

        - sum((h,b), mu(h,b) * dem_b(h,b) )$sw_elas
;

eq_dembin_lim(b,h)$sw_elas..
        qlim_bin(h,b) =g= dem_b(h,b) ; 

eq_capcon(s,f,pid,h)$genfeas(s,f,pid,h).. 
        capfac(f) * plantdata(s,f,pid,"cap") =g= X(s,f,pid,h);

eq_demcon(h).. sum((s,f,pid)$genfeas(s,f,pid,h),X(s,f,pid,h)) 
                =g= d(h)$(not sw_elas) + sum(b,dem_b(h,b))$sw_elas;

eq_carboncap(s)$SwCap(s).. 
        (1-psi(s)) * ebar_cap(s) 
        =g= 
        sum((f,pid,h)$genfeas(s,f,pid,h),
                plantdata(s,f,pid,"hr") * emit(f) * X(s,f,pid,h)
                );

eq_ratestandard(s)$SwTPS(s).. 
                (1-psi(s)) * ebar_rate(s) 
                * sum((f,pid,h)$genfeas(s,f,pid,h),X(s,f,pid,h)) 
                =g= 
                sum((f,pid,h)$genfeas(s,f,pid,h),
                        plantdata(s,f,pid,"hr") * emit(f) * X(s,f,pid,h)
                   );
;

model ele_dispatch /all/;

solve ele_dispatch using LP minimizing Z;

execute_unload 'primal.gdx' ; 

equation obj_dual, zpc_dual ; 

variable z_dual ; 

positive variables
lambda_k, lambda_d ; 

obj_dual.. z_dual =e=
        sum((h),lambda_d(h) * d(h))
        - sum((s,f,pid,h)$genfeas(s,f,pid,h), 
                                        lambda_k(s,f,pid,h) * capfac(f) * plantdata(s,f,pid,"cap") )
;

zpc_dual(s,f,pid,h)$genfeas(s,f,pid,h)..
             (plantdata(s,f,pid,"onm") + plantdata(s,f,pid,"hr") * fcost(f))
           +  lambda_k(s,f,pid,h) =g= lambda_d(h) ; 

model dual /obj_dual, zpc_dual/ ; 

solve dual using LP maximizing z_dual ; 

execute_unload 'dual.gdx' ; 

model mcp_ele 
/
eq_demcon.lambda_d,
eq_capcon.lambda_k,
zpc_dual.x
/ ; 

solve mcp_ele using mcp ; 


$exit


pbar(h) = eq_demcon.m(h) ; 
mu(h,b) = pbar(h) * (sum(bb$[bb.val<=b.val], qlim_bin(h,bb) ) / qbar(h) ) ** elasticity; 

sw_elas = 1 ; 

solve ele_dispatch using LP minimizing Z;

execute_unload 'eledat.gdx';




