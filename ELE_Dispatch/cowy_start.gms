$Title Simple Hourly Dispatch LP
* Maxwell Brown

$if not set sw_season $setglobal sw_season "day"
$if not set sw_co_reduction $setglobal sw_co_reduction 10
$if not set sw_wy_reduction $setglobal sw_wy_reduction 7.5


*Policy Switches
Scalar CAPTrade "Switch to turn on [1] or off [0] a carbon cap with trading" /0/;
Scalar TPSTrade "Switch to turn on [1] or off [0] a carbon TPS with trading" /0/;

*=======================
*  Begin sets
*=======================

set h "hour" /h1*h24/;
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

