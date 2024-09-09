
set i /coal, gas, wind/,  
    r /r1, r2/,
    h /1*24/ ; 

* have normally had simple examples like this...
parameter c(i)
/
Coal	20
Gas	15
Wind	1
/ ;


* can also store data in a .gms file:
parameter cf(h) "wind capacity factor"
* include an already-formatted-for-gams table
$include wind.gms
;

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

