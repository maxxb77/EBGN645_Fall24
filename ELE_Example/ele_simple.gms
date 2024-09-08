
set i /coal, gas, wind/,  
    r /r1, r2/,
    h /1*24/ ; 

* have normally made simple examples like this...
parameter c(i)
/
Coal	20
Gas	15
Wind	1
/ ;


* can also store data in a .gms file:
parameter cf(h) "wind capacity factor"
$include wind.gms
;

* note use of table in this next block
