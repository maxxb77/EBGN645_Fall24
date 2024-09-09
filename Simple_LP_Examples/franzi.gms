$TITLE Franzi 

*-- sets/indices --

set t "time" /0*30/
*    t0(t) "our first year" /0/
    t0(t) "our first year"
    tlast(t) "our last year"
 ; 
 
 alias(t,tt) ; 

set tprev(t,tt) "tt is the previous year to t" ; 

tprev(t,tt)$[tt.val = t.val-1] = yes ; 

t0(t)$[t.val=smin(tt,tt.val)] = yes ; 
tlast(t)$[t.val=smax(tt,tt.val)] = yes ; 

*--- Data ---
scalar 
    p  "sale price ($/tree)" /100/ 
    g  "annual growth rate for trees (%)" /1.021/
    s0 "initial stock of trees (trees)" /100/ 
    r  "discount rate (%)" /0.02/ 
; 

positive variable 
    H(t) "amount harvested in each year (trees)"
    S(t) "stock of trees available in each (trees)"
;

variable Z "objective function value" ; 

equation 
    objfn "target of our optimization"  
    eq_harvestlimit(t) "cant harvest more than what is in the stand"
    eq_stock_tracking(t) "inter-year tracking of available tree stock"
;

* three options for constraints:
* =e= : equal to
* =l= : less than
* =g= : greater than 
* net present of profit is equal to net profit * discount rate * price * harvest
* 1 / ((1+r)**t.val) ; 
scalar sw_salvage "enable or disable salvage value" /0/ ; 
objfn.. Z =e= sum(t, (1 / ((1+r)**t.val)) * p * H(t))
              + sum(t$tlast(t), sw_salvage * S(t) ); 

eq_harvestlimit(t)..
    S(t) =g= H(t)
;

eq_stock_tracking(t)..
    S(t) 
    =e= 
    
    s0$t0(t)

    + (g * sum(tt$tprev(t,tt), S(tt) - H(tt) ))$[not t0(t)]

    - H(t)$[tlast(t)$sw_salvage]
    
;

model franzi /all/;

solve franzi using lp maximizing z ; 

execute_unload 'franzi_data.gdx' ; 

$exit

parameter rep "reporting parameter";

set disc_loop /0,5,10,15,20,25,30/ ; 
set growth_loop /0,1,5,10,15,20/ ; 

alias(disc_loop,dl),(growth_loop,gl) ; 

loop(dl,
    loop(gl,
    r = dl.val / 100; 
    g = gl.val / 100; 
    solve franzi using lp maximizing z ; 
    rep(dl,gl,"H",t) = H.l(t) ; 
    );
); 


execute_unload 'franzi_data.gdx' ; 