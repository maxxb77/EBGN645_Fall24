$title Finite-Horizon Ramsey Growth Model as an NLP

SET    t   time period /0*50/;

* The model logic is different in the first and last periods of the
* model, so we define some subsets here to keep track of these
* periods.

SETS
        tfirst(t)       first period
        tlast(t)        last period
;

tfirst(t) = yes$(ord(t) eq 1);
tlast(t) = yes$(ord(t) eq card(t));

PARAMETERS
*  assumed...
        g       labor growth rate                       /0.023/
        delta   capital depreciation                    /0.04/
        i0      base year investment                    /0.3/
        c0      base year consumption                   /0.27/
        b       value share of capital in C-D prod fn   /0.65/
        eta     inverse intertemporal elasticity        /2/

*calibrated parameters...
        k0      initial capital stock
        l0      initial labor supply
        l(t)    labor supply
        a       prod fn scale parameter
        rho     pure time preference rate
        qref(t) steady-state output index
        beta(t) discount factor
;