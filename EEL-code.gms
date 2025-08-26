$title Epi-Eco LEWIE Model
*by Edward Whitney

Option limrow=0 ;
Option limcol=0 ;
Option DECIMALS=5 ;

*==============================================================================
*===================== Choose Shock and Trade Scenarios =======================
*==============================================================================

*Define Shock parameters (equal to zero by default)
parameter MDAshock equal to 1 if simulating MDA shock  / 0 /
          FMPshock equal to 1 if simulating FMP shock  / 0 /
          TFPshock equal to 1 if simulating TFP shock  / 0 /
;

*Choose which shock(s) to activate by setting the corresponding shock parameter(s) equal to 1
*MDAshock = 1 ;
*FMPshock = 1 ;
*TFPshock = 1 ;

*Define Trade parameter (equal to 1 for preferred specification)
parameter RoWtrade equal to 1 if Fish is traded RoW and zero if traded Local / 1 /
;

*Is Fish traded Locally? If so, activate next line
*RoWtrade = 0 ;

*==============================================================================
*================== Initial Set-up ============================================
*==============================================================================

*define t in terms of years
set      t   time step representing one year    / 1*10 /;

* Name the sets that will be used:
sets
ac all accounts
g(ac) goods
f(ac) factors
h(ac) households
var  variable names
;

*define a parameter that contains all the data from the input sheet:
Parameters Alldata(*,*,*,*,h);

* define model tolerance threshold (larger than this value and the equation is deemed infeasible)
$setglobal model_tolerance "1.0E-10"

*==============================================================================
*======================== Read in data input sheet ============================
*==============================================================================

* call in the data file
$gdxin EEL-data.gdx
$load AC G F H VAR ALLDATA
display ac, g, f, h, var, ALLDATA ;

* This option controls the decimals and display format
option alldata:4:4:1;
display alldata;

* the phantom element "null" can be put in a set to avoid leaving the
* set empty in some simulations (GAMS can't handle empty sets)
$phantom null

* Define set aliases
*======================================================
alias (ac, aca) ;
alias (g,gg,ggg,gfac)
      (h,hh)
      (f,fa) ;

*======================================================================
*================== DEFINE INITIAL PARAMETERS =========================
*======================================================================

* Get raw values for input variables and some parameter values
parameter
* expenditure function parameters:
     xleshare(g,h)         expenditure share of good g by household h
     xleshare_se(g,h)      expenditure share standard error
     xlemin(g,h)           minimum expenditure on g

* production function parameters:
     xlidsh(gg,g,h)        intermediate demand share for of gg for production of g by h
     xlqp(g,h)             quantity of g produced by household h
     xlqp_se(g,h)          standard error for xlqp
     xlfshare(g,f,h)       exponent on factor f in produciton of g
     xlfshare_se(g,f,h)    standard error on beta (gfh)
     xlpshift(g,h)         shift parameter on production of g
     xlpshift_se(g,h)      standard error on acobb

* other parameters:
     xlendow(f,h)          endowment of factors in the economy
     xlROCendow(f,h)       endowment of factors outside the economy
     xlROWendow(f,h)       endowment of factors outside the country
     xlTRINsh(h)           cash transfers given to other households (share of income)
     xlTROUTsh(h)          cash transfers received from other households  (share of expenditures)
     xlTRINsh_se(h)         standard error of cash transfers given to other households (share of income)
     xlTROUTsh_se(h)        standard error of cash transfers received from other households  (share of expenditures)
     xlsavsh(h)         share of income going to informal savings
     xlsavsh_se(h)       standard error of share of income going to informal savings
     xllabexp(h)           not sure what this is and why there's a negative value
     xlexpoutsh(h)         share of expenditures on outside goods
     xlremit(h)            level of remittances
     xlothertransfers(h)   level of exogenous transfers
     xlnhh(h)              number of households represented by this
     xlhhinc(h)            mean household income
     xlhhexp(h)            mean household expenditures
     xlhhsize(h)           mean household size
     xlrevsh_vil(g,h)        share of business in the village
     xlrevsh_zoi(g,h)        share of business in the zoi
     xlrevsh_rol(g,h)        share of business in the rest of lesotho
     xlrevsh_row(g,h)        share of business in the row
     xlVA2IDsh(g,gg,h)     for each dollar of VA how much ID was consumed
     xlwrkagepop(h)           number of working-age ppl in the household
;

* expenditure parameters
xleshare(g,h) = alldata("eshare",g,"","",h) ;
xleshare_se(g,h) = alldata("eshare_se",g,"","",h) ;
xlemin(g,h) = alldata("emin",g,"","",h) ;

* production parameters
xlqp(g,h)      = alldata("qp",g,"","",h) ;
xlqp_se(g,h)   = alldata("qp_se",g,"","",h) ;
xlidsh(g,gg,h) = alldata("idsh",g,gg,"",h) ;
xlfshare(g,f,h) = alldata("fshare",g,"",f,h) ;
xlfshare_se(g,f,h) = alldata("fshare_se",g,"",f,h) ;
xlpshift(g,h)= alldata("pshift",g,"","",h) ;

* transfers and savings
xlTROUTsh(h) = alldata("transfoutsh","","","",h) ;
xlTRINsh(h) = alldata("transfinsh","","","",h) ;
xlTROUTsh_se(h) = alldata("transfoutsh_se","","","",h) ;
xlTRINsh_se(h) = alldata("transfinsh_se","","","",h) ;
xlsavsh(h) = alldata("savsh","","","",h) ;
xlSAVsh_se(h) = alldata("savsh_se","","","",h) ;
xlendow(f,h) = alldata("endow","","",f,h) + alldata("zoiendow","","",f,h) ;
xlROCendow(f,h) = alldata("ROCendow","","",f,h) ;
xlROWendow(f,h) = alldata("ROWendow","","",f,h) ;
xlexpoutsh(h) = alldata("exproles","","","",h) ;
xlremit(h)  =  alldata("remits","","","",h) ;
xlothertransfers(h)  =  alldata("NONSCtransfers","","","",h) ;
xlnhh(h) = alldata("HHNum","","","",h) ;
xlhhinc(h) = alldata("HHinc","","","",h) ;
xlhhexp(h) = alldata("HHexp","","","",h) ;
xlhhsize(h) = alldata("HHsize","","","",h) ;
xlrevsh_vil(g,h) = alldata("revsh_vil",g,"","",h) ;
xlrevsh_zoi(g,h) = alldata("revsh_zoi",g,"","",h) ;
xlrevsh_rol(g,h) = alldata("revsh_rol",g,"","",h) ;
xlrevsh_row(g,h) = alldata("revsh_row",g,"","",h) ;
xlwrkagepop(h)      = alldata("wrkagepop","","","",h)
;

*==============================================================================
*============================= Define Most Sets ===============================
*==============================================================================

* subsets and aliases
*=====================
set
* factor subsets
     fk(f)     fixed factors /CAPITAL,  LAND /
     ft(f)     tradable inputs / LABOR, INPUT /
     ftv(f)    factors tradable in the village /LABOR,  null /
     ftw(f)    factors tradables in the rest of the world / INPUT /
     fpurch(f) purchased factors /INPUT /

* goods subsets
     gtv(g)    goods tradable in the village / ser, ret, meat, crop /
     gtw(g)    goods tradable with the rest of the world / out, palmoil, fish  /

     gm(g)     village-traded goods with imperfect substitutes / fish /
     gd(g)     village-traded goods withOUT imperfect substitutes / ser, ret, meat, crop /

     gp(g)     goods that are produced by at least one hh group
                  /  palmoil, ser, ret, meat, fish, crop /
     gag(g)    ag goods / palmoil, crop, meat, fish /
     gnag(g)   non ag goods / ser, ret, null /
*Make a subset of good for those goods that have a stock.
     gstk(g)   produced goods that have a stock (e.g. the fish good has a stock level)   / fish /
;

* accounts not in the matrix
sets
     v        villages / HHs  /

     maphv mapping household to their village / (FishNpoor , FishPoor, NfishNpoor, NfishPoor).HHs
/
;

*==============================================================================
*================= REVISE SETS AS NEEDED BASED ON =============================
*================= MARKET STRUCTURE FOR SIMULATION ============================
*==============================================================================

* Include "fish" in locally traded goods for Local scenario
gtv("fish")$(RoWtrade = 0) = yes ;

* Exclude "fish" in RoW traded goods for Local scenario
gtw("fish")$(RoWtrade = 0) = no ;

*==============================================================================
*============================= ELASTICITY =====================================
*==============================================================================

* set the elasticity of supply of labor
$setglobal hlse 100

*TRADE ELASTICITY
*Define and set trade elasticity for Armingtom function(s):
parameter trade_elas(g)  trade elasticity for CES aggregation function goods with imperfect substitutes ;
trade_elas(g) = 1;
trade_elas("fish")      = 8;

*CES UTILITY ELASTICITIES:
*Define and set elasticities of substitutions for utility functions:
parameter good_elas(h)    elasticity of substitution for goods for HHs ;

*Define the elasticity of substitution for goods in each HH's utility function
good_elas(h)     = 3 ;


*==============================================================================
*======================== INITIAL FISHING PARAMETERS ==========================
*==============================================================================

*SEARCH COST PARAMETER (n_sc):
parameters       n_sc     search cost parameter (exponent in denominator)

*GROWTH RATE  (for stock difference equation):
                 gr(gstk)              growth rate

                 stock_dr(g)        stocks for the production functions that have stocks (set =1 for those don't)
                 stockbeta(g)    stocks production coefficients for the production functions that have stocks (set =0 for those don't)

*Biological information (growth rate is defined above):
                 k_dr(gstk)                  carrying capacity

                 x_bl(g)                    baseline stock level (for calculating pc in results output)
                 harvestkg_bl(gstk)                       harvest at baseline of stock good gstk
                 harvestkg_ts(gstk,t)                       harvest at time t of stock good gstk
                 harvestkg_d(gstk,t)                     difference in harvest at time t (compared to baseline)
                 harvestkg_pc(gstk,t)                     percent change in harvest at time t (compared to baseline)
                 idsh_fish_record_ts(g,gg,h,t)       this records the idsh for fishing sector (to see if it changes with stock changes)
                 idsh_fish_record_pc(gstk,gg,h,t)    percentage changes of idsh for fishing sector
                 search_cost_record_ts(g,t)          records search cost coeff. used to adjust PVA in fishing sector
                 PCdomestic_import_ratio_ts(g,t)     ratio of domestic production to imports

                 tt                                      counter for fish transition equation
;

*ASSIGN INITIAL VALUES
n_sc             = 2;

*growth rate (annual)
gr("fish")      = 1.06;

display gr ;

*stock betas
stockbeta(g)$gp(g)  = 1;
stockbeta("fish")  = 0.645;

*Stock size parameters
stock_dr(g)$gp(g)    =  1;

*==============================================================================
*=========================== DEFINING THE MODELS ==============================
*==============================================================================

* MODEL STARTS HERE
* ======================================================
* Variables and parameters for the LEWIE component
* ======================================================

nonnegative variables
* production
     QP(g,h)        quantity produced of a good by a household
     FD(g,f,h)      factor demand of f in production of g
     ID(g,gg,h)     intermediate demand of gg for production of g
     QVA(g,h)       quantity of value added created in the production process

     HFD(f,h)       factor demand in the household
     HFSUP(f,h)     Effective labor supply from the household (elastic endowment)
     LABTIME_BLvar(f,h) Labor time supplied by the household
     VFD(f)       initial factor demand in the village

     R(g,f,h)       rent for fixed factors
     WV(f)        wage at the village level

* consumption
     QC(g,h)        quantity consumed of g by h
     Y(h)           nominal household income
     RY(h)          real household income
     CPI(h)         consumer price index

* prices
     P(g)           price of good
     PVA(g,h)       price of value added net of intermediate inputs as seen by the household


* transfers
     TRIN(h)        tranfered in - received by a household
     TROUT(h)       transfers out - given by a household
     SAV(h)         household savings
     EXPROC(h)      household expenditures out of the zoi

*Armington varialbes:
     QP_COMPOSITE(g)     Quantity of composite good produced (using domestic goods and imports)
     P_COMPOSITE(g)      Price of composite good
     IMPORTS(g)          Level of imports for goods that have both domestic production and imports
;


variables
* trade
     HMS(g,h)  household marketed surplus of good g
     VMS(g)    village marketed surplus of good g

     HFMS(f,h) factor marketed surplus from the household
     VFMS(f)   factor marketed surplus out of the village
;

parameters
*Production - Cobb-douglas
     pshift(g,h)         production shift parameter for the CD
     fshare(g,f,h)       factor share parameter for the CD
     vash(g,h)           share of value added
     idsh(g,gfac,h)      intermediate input share
     tidsh(g,h)          total intermediate input share (1-vash)
     sumbetas_m(g,h)     sum of betas for model solve (equal to value in each draw)

*Consumption
     eshare(g,h)         expenditure share parameters in the LES
     emin(g,h)           minimal expenditure in the LES
     exinc(h)            exogenous income of household
     vmsfix(g)           fixed marketed surplus at the village level
     utility_sh_m(g,h)    utility share (equal to value for each draw)
     delta_m(g)          share parameter for the armington CES (equal to value for each draw)
     armg_shifter_m(g)   shift coefficient for the armington CES (equal to value for each draw)
     rho(g)                   CES coefficient for the armington CES aggregation function
     p_import_m(g)           the price of imports of good type g (equal to value for each draw)

* factor endowments for fixed factors
     fixfac(g,f,h)       fixed factors
     vfmsfix(f)          factors fixed at the Village level (family labor)
     endow(f,h)          endowment of factors
     hfsupzero(f,h)      initial labor supply

* Factor supply
     hfsupel(f,h)        factor supply elasticity from household
     labtime_bl(f,h)     Baseline labor time supplied by the household (PARAMETER IS FOR SIMS MODEL ONLY... assigned values from BL Solve)

* Transfers
     troutsh(h)          share of transfers in the households expenditures
     exprocsh(h)         share of expenditures outside of the zoi
     savsh(h)            share of income saved
     trinsh(h)           share of total transfers received by a given household

* Fish model parameters for model solve
     stock(g)      stock parameter for model solve (equal to value in each draw)
     theta_fish1(g,f,h)    theta parameter for model solve (equal to value in each draw)
     sumbetas(g,h)       Sum of production betas (inclusive of stock beta)

;

* ==================================================================
* Now parameters for the Epi component and baseline model variables
* ==================================================================

NONNEGATIVE VARIABLES
fishlabtime                  Total time supplied to fishing sector (solved in baseline SS model)
fishlabtime_hh(h)             Household time supplied to fishing sector (solved in baseline SS model)
Infectrate_hh(h)              Infection Rate for household (solved in baseline SS model)
extime_blvar(h)               Baseline exposure time for the household
Ystate                        Y state variable (solved in baseline SS model)
beta                          Baseline value for beta
gamma                         Baseline value for gamma
;

parameters
         mu              mortality rate of snails
         bmin            Minimum value for beta parameter
         bmax            Maximum value for beta parameter
         gmin            Minimum value for gamma parameter
         gmax            Maximum value for gamma parameter
         ymed            Median value for yt fraction
         yslope          Slope value for yt fraction
         alpha           The portion of labor productivity lost due to infection
         vareps_hh(h)    Household E-C risk parameter
         sigma_scalar    Sigma Scalar
         lnagg_m         Logged Aggregate Output Per Capita (equal to value for each draw)
         lnagg_med       Median value for logged aggregate output pc
         lnagg_s         Slope value for logged aggregate output pc
         blextime(h)     Baseline exposure time for the household (PARAMETER IS FOR SIMS MODEL ONLY... assigned values from BL Solve)

* Infection rate for SIMS model
      Infectrate_hh_sims(h)    Infection rate for household in SIMS model

        shhh(h) household share of total household population
;


* Treatment that increases mortality of parasite in human host
*===============================================================

parameter
         mdaeffcov             model value for mdaeffcov(t)
        mdaeffcov_ts(t)       Combined efficacy and converage for MDA
;

*==============================================================================
*=== Now Define Variables and Parameters for finding Quasi-SS values over time
*==============================================================================

NONNEGATIVE VARIABLES
Infectrate_hh_ts_t(h)              Temp Infectrate_hh_ts
Ystate_ts_t                        Temp Y state variable
;

parameters       beta_ts_t beta in time t (solved for in loop in main gms file)
                 gamma_ts_t gamma in time t (solved for in loop in main gms file)
                 extime_ratio_t(h) extime_ratio at time t (solved for in loop in main gms file)
;

***LEWIE Component Equations
Equations
* prices
     EQ_PVA(g,h)         privet value added equation

* production
     EQ_FDCOBB(g,f,h)    factor demands cobb douglas
     EQ_FDPURCH(g,f,h)   factor demands for purchased inputs - constrained or not
     EQ_QVACOBB(g,h)     quantity VA produced cobb douglas
     EQ_QP(g,h)          quantity produced from QVA and ID
     EQ_ID(gg,g,h)       quantity of g ID needed for QP of gg

* consumption
     EQ_QC(g,h)          quantity consumed

* income
     EQ_Y(h)          full income constraint for the household
     EQ_CPI(h)           consumer price index equation
     EQ_RY(h)            real household income equation

* transfers
     EQ_TRIN(h)          inter household transfers in (received)
     EQ_TROUT(h)         interhousehold transfers out (given)

* exogenous expenditures
     EQ_SAV(h)           savings (exogenous rate)
     EQ_EXPROC(h)        expenditures outside of the zoi (exogenous rate)

* goods market clearing
     EQ_HMKT(g,h)        qty clearing in each household
     EQ_VMKT(g)        market clearing in the village
     EQ_VMKTfix(g)       price definition in the village

* factor market clearing
     EQ_HFD(f,h)         total household demand for a given factor
     EQ_FCSTR(g,f,h)     fixed factors constraint

     EQ_HFSUP_BL(f,h)    household elastic supply (For BL ONLY)
     EQ_LABTIME_BL(f,h)  household supply of labor time at baseline
     EQ_HFSUP_SIMS(f,h)  household elastic supply (For SIMS ONLY)

     EQ_HFMKT(f,h)       tradable factor clearing in the household
     EQ_VFMKT(f)         tradable factors clearing in the village
     EQ_VFMKTfix(f)       wage determination for tradable factors clearing in the village

*Armington function equations
     EQ_QP_COMPOSITE(gm)     defines the composite commodity for goods with imperfect substitute imports
     EQ_P_COMPOSITE(gm)      defines the price of the composite commodity
     EQ_IMPORTS(gm)          defines the level of imports
;

***Epi Component Equations Baseline Model
EQUATIONS
EQ_Infss(h)     Equation for Steady-state value of Infection Rate for household
EQ_Yss          Equation for Steady-State value for Y state var
EQ_Lfishhh(h)   Equation for Steady-State Household Fishing Labor time
EQ_Lfish        Equation for Steady-State Total Fishing Labor time
EQ_beta         Equation for SS beta
EQ_gamma         Equation for SS gamma
;

***Epi Component Equations Sims Model
EQUATIONS
EQ_Lfishhh_ts(h)   Equation for Steady-State Household Fishing Labor time
EQ_Lfish_ts        Equation for Steady-State Total Fishing Labor time
EQ_beta_ts         Equation for SS beta
EQ_gamma_ts         Equation for SS gamma
;

***Epi Component Equations For Quasi-SS Model
EQUATIONS
EQ_Infss_ts(h)     Equation for Quasi-SS value of HH Infection Rate
EQ_Yss_ts          Equation for Quasi-SS value for Y state var
;

*================================================================================
*==================== LEWIE COMPONENT EQUATIONS =================================
*================================================================================

EQ_PVA(g,h)..
     PVA(g,h) =E= P(g)- sum(gg,idsh(g,gg,h)*( P(gg)$(not gm(gg)) + P_COMPOSITE(gg)$gm(gg)) )
;

EQ_QVACOBB(g,h)..
     QVA(g,h) =E= pshift(g,h)*prod(f,FD(g,f,h)**(fshare(g,f,h)))*stock(g)**stockbeta(g)
;

EQ_FDCOBB(g,f,h)$gp(g)..
     FD(g,f,h)*(R(g,f,h)$fk(f) + WV(f)$ftw(f) + WV(f)$ftv(f))*sumbetas(g,h)
      =E= PVA(g,h)*QP(g,h)*[fshare(g,f,h) + theta_fish1(g,f,h)*stockbeta(g)]
;

EQ_QP(g,h)$vash(g,h)..
     QP(g,h) =E= QVA(g,h)/vash(g,h)
;

EQ_ID(g,gfac,h)..
     ID(g,gfac,h) =E= QP(g,h)*idsh(g,gfac,h)
;

* CONSUMPTION AND INCOME
EQ_QC(g,h)$(utility_sh_m(g,h)>0)..
      QC(g,h) =E= [ ( P(g)$(not gm(g)) + P_COMPOSITE(g)$gm(g) )**(1-good_elas(h)) *  utility_sh_m(g,h)**(good_elas(h)-1) * (Y(h)-TROUT(h)-SAV(h)-EXPROC(h))] /
                  [ ( P(g)$(not gm(g)) + P_COMPOSITE(g)$gm(g) ) * sum(gg, ( ( P(gg)$(not gm(gg)) + P_COMPOSITE(gg)$gm(gg) )**(1-good_elas(h)) *  utility_sh_m(gg,h)**(good_elas(h)-1))) ] ;

;

EQ_Y(h)..
     Y(h) =E= sum((g,fk),R(g,fk,h)*FD(g,fk,h))
            + sum(ftv, WV(ftv)*HFSUP(ftv,h))
            + sum(ftw, WV(ftw)*HFSUP(ftw,h))
            + exinc(h)
;

EQ_CPI(h)..
     CPI(h) =e= sum(g, ( P(g)$(not gm(g)) + P_COMPOSITE(g)$gm(g) ) *[ QC(g,h)/sum(gg,QC(gg,h)) ] )
;

EQ_RY(h)..
     RY(h) =e= Y(h) / CPI(h)
;

* Transfers given away
EQ_TROUT(h)..
     TROUT(h) =E= troutsh(h)*Y(h) ;
;

EQ_SAV(h)..
     SAV(h) =E= savsh(h)*Y(h) ;
;

EQ_EXPROC(h)..
     EXPROC(h) =E= exprocsh(h)*Y(h) ;
;

* MARKET CLEARING FOR GOODS
EQ_HMKT(g,h)..
      HMS(g,h) =E=  QP(g,h)$(not gm(g)) - QC(g,h) - sum(gg,ID(gg,g,h))
;

EQ_VMKT(g)..
      VMS(g) =E= QP_COMPOSITE(g)$gm(g)+ sum(h, HMS(g,h))
;

EQ_VMKTfix(gtv)..
     VMS(gtv) =E= vmsfix(gtv)
;

* FACTOR MARKET CLEARING
EQ_HFD(f,h)..
     HFD(f,h) =e= sum(g, FD(g,f,h))
;

EQ_FCSTR(g,fk,h)..
     FD(g,fk,h) =E= fixfac(g,fk,h)
;

EQ_HFMKT(ft,h)..
     HFMS(ft,h) =E= HFSUP(ft,h) - sum(g, FD(g,ft,h))
;

EQ_VFMKT(ft)..
     VFMS(ft) =E= sum(h, HFMS(ft,h))
;

* Elastic Labor supply for BASELINE
EQ_HFSUP_BL(ft,h)..
     HFSUP(ft,h)$(not hfsupzero(ft,h))
        +
     (HFSUP(ft,h)/(labtime_blvar(ft,h)*(1-Infectrate_hh(h)*alpha))  - [sum(v$maphv(h,v),WV(ft)**hfsupel(ft,h))$ftv(ft) + (WV(ft)**hfsupel(ft,h))$(ftw(ft))] )$hfsupzero(ft,h)
     =e= 0
;

*** Since Labtime_bl(ft,h) is unknown
EQ_LABTIME_BL(ft,h)..
     LABTIME_BLvar(ft,h)$(not hfsupzero(ft,h))
        +
     [LABTIME_BLvar(ft,h)*(1-Infectrate_hh(h)*alpha) - hfsupzero(ft,h)]$hfsupzero(ft,h)
     =e= 0
;

* Elastic Labor supply for SIMULATIONS
EQ_HFSUP_SIMS(ft,h)..
     HFSUP(ft,h)$(not hfsupzero(ft,h))
        +
     (HFSUP(ft,h)/(labtime_bl(ft,h)*(1-Infectrate_hh_sims(h)*alpha))  - [sum(v$maphv(h,v),WV(ft)**hfsupel(ft,h))$ftv(ft) + (WV(ft)**hfsupel(ft,h))$(ftw(ft))] )$hfsupzero(ft,h)
     =e= 0
;

* Equations related to the Armington function:
     EQ_QP_COMPOSITE(gm)..
          QP_COMPOSITE(gm) =e=  armg_shifter_m(gm)* ( delta_m(gm)*sum(h,QP(gm,h))**(rho(gm))
                               + (1-delta_m(gm))*IMPORTS(gm)**(rho(gm)))**(1/rho(gm));
     EQ_IMPORTS(gm)..
          IMPORTS(gm) =e= sum(h,QP(gm,h))*(  (p_import_m(gm)/P(gm))*(delta_m(gm)/(1-delta_m(gm))  )  )**(1/(rho(gm)-1));

     EQ_P_COMPOSITE(gm)..
          P_COMPOSITE(gm) =e= ( p(gm)*sum(h,QP(gm,h)) + p_import_m(gm)*IMPORTS(gm) )  / QP_COMPOSITE(gm);

* FACTOR WAGE DETERMINATION
EQ_VFMKTFIX(ftv)..
     VFMS(ftv) =E= vfmsfix(ftv)
;

*==============================================================================
*=============== EPI COMPONENT EQUATIONS: Baseline Model ======================
*===== (Infectrate_hh(h) is a variable to be solved for) ======================
*==============================================================================


EQ_infss(h)..
gamma*Infectrate_hh(h) =e= beta*vareps_hh(h)*Ystate*(1-Infectrate_hh(h)) ;

EQ_Yss..
Ystate*mu =e= (1-Ystate)*beta*sum(h,shhh(h)*vareps_hh(h)*Infectrate_hh(h)) ;

EQ_Lfishhh(h)..
fishlabtime_hh(h) =e= FD("fish","labor",h)/(1-Infectrate_hh(h)*alpha) ;

EQ_Lfish..
fishlabtime =e= sum(h,fishlabtime_hh(h)) ;

EQ_beta..
beta =e= (1 + exp((log(sum(g,sum(h,QP(g,h)))*1000000/sum(h,xlnhh(h)*xlhhsize(h)))-lnagg_med)/lnagg_s))**(-1)*(bmax-bmin) + bmin ;

EQ_gamma..
gamma =e= (1 + exp(-(log(sum(g,sum(h,QP(g,h)))*1000000/sum(h,xlnhh(h)*xlhhsize(h)))-lnagg_med)/lnagg_s))**(-1)*(gmax-gmin) + gmin ;


*==============================================================================
*=============== EPI COMPONENT EQUATIONS: Sims Model ==========================
*= (Infectrate_hh_sims(h) is a parameter, value identified in QSSA model) =====
*==============================================================================

EQ_Lfishhh_ts(h)..
fishlabtime_hh(h) =e= FD("fish","labor",h)/(1-Infectrate_hh_sims(h)*alpha) ;

EQ_Lfish_ts..
fishlabtime =e= sum(h,fishlabtime_hh(h)) ;

EQ_beta_ts..
beta =e= (1 + exp((log(sum(g,sum(h,QP(g,h)))*1000000/sum(h,xlnhh(h)*xlhhsize(h)))-lnagg_med)/lnagg_s))**(-1)*(bmax-bmin) + bmin ;

EQ_gamma_ts..
gamma =e= (1 + exp(-(log(sum(g,sum(h,QP(g,h)))*1000000/sum(h,xlnhh(h)*xlhhsize(h)))-lnagg_med)/lnagg_s))**(-1)*(gmax-gmin) + gmin ;

*==============================================================================
*=============== EPI COMPONENT EQUATIONS: QSSA MODEL ==========================
*==============================================================================

EQ_infss_ts(h)..
(gamma_ts_t + mdaeffcov)*Infectrate_hh_ts_t(h) =e= beta_ts_t*extime_ratio_t(h)*vareps_hh(h)*Ystate_ts_t*(1-Infectrate_hh_ts_t(h)) ;

EQ_Yss_ts..
Ystate_ts_t*mu =e= (1-Ystate_ts_t)*beta_ts_t*sum(h,shhh(h)*extime_ratio_t(h)*vareps_hh(h)*Infectrate_hh_ts_t(h)) ;

Infectrate_hh_ts_t.lo(h) = 0 ;
Infectrate_hh_ts_t.up(h) = 1 ;
Ystate_ts_t.lo = 0 ;
Ystate_ts_t.up = 1 ;


*==========================================================
*================== MODEL STATEMENTS ======================
*==========================================================

model BL_MODEL Baseline Epi-LEWIE model /
EQ_PVA.PVA
EQ_QVACOBB.QVA
EQ_FDCOBB.FD
EQ_QP.QP
EQ_ID.ID
EQ_QC.QC
EQ_Y.Y
EQ_HMKT.HMS
EQ_VMKT.VMS
EQ_VMKTfix.P
EQ_HFD.HFD
EQ_FCSTR.R
EQ_HFMKT.HFMS
EQ_VFMKT.VFMS
EQ_VFMKTfix.WV
EQ_TROUT.TROUT
EQ_SAV.SAV
EQ_EXPROC.EXPROC
EQ_QP_COMPOSITE.QP_COMPOSITE
EQ_IMPORTS.IMPORTS
EQ_P_COMPOSITE.P_COMPOSITE

* elastic factor supply from the household
EQ_HFSUP_BL.HFSUP
EQ_LABTIME_BL.LABTIME_BLvar

EQ_CPI.CPI
EQ_RY.RY

EQ_infss.Infectrate_hh
EQ_Yss.Ystate
EQ_Lfishhh.fishlabtime_hh
EQ_Lfish.fishlabtime
EQ_beta.beta
EQ_gamma.gamma
/;

* define tolerance for model solve
BL_MODEL.tolInfRep = %model_tolerance% ;


model SIMS_MODEL LEWIE Model for Simulations /
EQ_PVA.PVA
EQ_QVACOBB.QVA
EQ_FDCOBB.FD
EQ_QP.QP
EQ_ID.ID
EQ_QC.QC
EQ_Y.Y
EQ_HMKT.HMS
EQ_VMKT.VMS
EQ_VMKTfix.P
EQ_HFD.HFD
EQ_FCSTR.R
EQ_HFMKT.HFMS
EQ_VFMKT.VFMS
EQ_VFMKTfix.WV
EQ_TROUT.TROUT
EQ_SAV.SAV
EQ_EXPROC.EXPROC
EQ_QP_COMPOSITE.QP_COMPOSITE
EQ_IMPORTS.IMPORTS
EQ_P_COMPOSITE.P_COMPOSITE

* elastic factor supply from the household
EQ_HFSUP_SIMS.HFSUP

EQ_CPI.CPI
EQ_RY.RY

EQ_Lfishhh_ts.fishlabtime_hh
EQ_Lfish_ts.fishlabtime
EQ_beta_ts.beta
EQ_gamma_ts.gamma
/;

* define tolerance for model solve
SIMS_MODEL.tolInfRep = %model_tolerance% ;


***Quasi-SS Epi Model
Model qssa Model for Quasi-SS values of infection dynamics /
EQ_infss_ts.Infectrate_hh_ts_t
EQ_Yss_ts.Ystate_ts_t
/ ;

* define tolerance for model solve
qssa.tolInfRep = %model_tolerance% ;


*==============================================================================
*================== DEFINE ESSENTIAL CALIBRATION PARAMETERS ===================
*==============================================================================


parameter
* model result stat
modstat model stat (1 means solved ok)
modstat_dr model stat for each draw (1 means solved ok)

* meta-parameters with parameter draws
fshare_t(g,f,h)  unscaled draw the cobb-douglas factor shares
eshare_t(g,h)    unscaled draw of expenditure shares
qp_t(g,h)        unscaled draw of qp


* drawn parameters:
p_dr(g)          price of good
pv_dr(g,v)       price at village level
pz_dr(g)         price at zoi level
ph_dr(g,h)       price as seen from household
pva_dr(g,h)      price of value added
qva_dr(g,h)      quantity of value added
qp_dr(g,h)       quantity produced
tqp_dr(g)        total qty produced in the zoi
imports_dr(g)    total imports of g
ttqp_dr         total output of the zoi
ttqp_percap_dr  per capita output of the zoi
lnagg_dr        Logged total output (per capita)
fd_dr(g,f,h)     factor demand
id_dr(g,gg,h)    intermediate demand
pshift_dr(g,h)   cobb-douglas production shifter
fshare_dr(g,f,h) cobb-douglas factor shares
r_dr(g,f,h)      rent for fixed factors
wv_dr(f)       village-wide wage for tradable factors
vash_dr(g,h)     value-added share
idsh_dr(gg,g,h)  intermediate demand share
tidsh_dr(gg,h)   total intermediate input share (1-vash)
fixfac_dr(g,f,h) fixed factor demand
unemp_dr(f,h)    unemployment
unempsh_dr(f,h)  hh share of total unemployment
vfmsfix_dr(f)  factors fixed at the Village level
vmsfix_dr(g)   goods fixed at the Village level
exinc_dr(h)      exogenous income
endow_dr(f,h)    endowment
qc_dr(g,h)       level of consumption
tqc_dr(g)        total qc
eshare_dr(g,h)    consumption shares
y_dr(h)          nominal hh income
cpi_dr(h)        consumer price index of hh
ry_dr(h)         real hh income
emin_dr(g,h)     incompressible demand
trin_dr(h)       transfers in - received
trout_dr(h)      transfers out - given
trinsh_dr(h)     share of all transfers in the eco going to h
troutsh_dr(h)    share of yousehold h's income being given as transfers
hfd_dr(f,h)      factor demand of household h for factor f
vfd_dr(f)      village demand for factor f
hms_dr(g,h)      household marketed surplus of good g
vms_dr(g)      village marketed surplus of good g
hfms_dr(f,h)     household factor marketed surplus
vfms_dr(f)     village factor marketed surplus
savsh_dr(h)      savings rate
sav_dr(h)        savings level
exprocsh_dr(h)   outside-of-zoi expenditures rate
exproc_dr(h)     outside-of-zoi expenditures level
expzoish_dr(h)   outside-of-zoi expenditures level
fixfac_t(g,f,h)          first try at changing the fixfac parameter (gets changed if negative)
negfixfac(g,f,h) check for negative values for each sim and draw
fixfacsim_dr_t(g,f,h) fixed factor after the simulation (temp - not corrected)
fixfacsim_dr(g,f,h) fixed factor after the simulation (corrected if needed)
negfixfacnum     check the number of negatives that needed to be corrected for each sim
harvest1(gstk)      harvest of stock good gstk
;

*===============================================================================
*================== DEFINE PARAMETERS FOR BASELINE VALUES TO ===================
*=================== BE USED TO COMPUTE CHANGES OVER TIME ======================
*===============================================================================

parameter        p_bl(g)                 price
                 pva_bl(g,h)             price of value added
                 qva_bl(g,h)             quantity of value added
                 qp_bl(g,h)              quantity produced by h
                 tqp_bl(g)               total quantity produced of g
                Qp_composite_bl(g)    Quntity of the composite good produced by armington function
                P_composite_bl(g)     price of composite good
                Imports_bl(g)         imports used in the creation of the composite good
                 ttqp_bl              total production value in whole economy
                 ttqp_percap_bl       per capita output of the zoi
                 hqp_bl(h)             total qty produced by a household
                 fd_bl(g,f,h)          factor demand
                 id_bl(g,gg,h)         intermediate demand
                 pshift_bl(g,h)        cobb-douglas shifter
                 fshare_bl(g,f,h)      cobb-douglas shares
                 r_bl(g,f,h)           rent for fixed factors
                 wv_bl(f)            village-wide wage for tradable factors
                 vash_bl(g,h)          value-added share
                 idsh_bl(gg,g,h)       intermediate demand share
                 tidsh_bl(gg,h)        total intermediate input share (1-vash)
                 fixfac_bl(g,f,h)      fixed factor demand
                 exinc_bl(h)           exogenous income
                 endow_bl(f,h)         endowment
                 qc_bl(g,h)            level of consumption
                 eshare_bl(g,h)        consumption shares
                 y_bl(h)               income of household
                 mry_bl                mean income in economy (weighted by nhh)
                 rytheil_bl            theil index for income
                 cpi_bl(h)             cpi
                 vqc_bl(g)           village consumption
                 vcpi_bl            village cpi
                 cri_bl(f)           rent weighted index
                 ry_bl(h)              real income
                 ty_bl                income total
                 try_bl               real income total
                 emin_bl(g,h)          incompressible demand
                 trin_bl(h)            transfers in - received
                 trout_bl(h)           transfers out - given
                 sav_bl(h)             savings
                 exproc_bl(h)          expenditure rest of country
                 trinsh_bl(h)          share of all transfers in the eco going to h
                 troutsh_bl(h)         share of household h income being given as transfers
                 hfd_bl(f,h)           factor demand of household h for factor f
                 vfd_bl(f)           village demand for factor f
                 zfd_bl(f)             zoi demand for factor f
                 hms_bl(g,h)           household marketed surplus of good g
                 vms_bl(g)           village marketed surplus of good g
                 hfms_bl(f,h)          household factor marketed surplus
                 vfms_bl(f)          village factor marketed surplus
                 vfmsfix_bl(f)       factors fixed at the Village level (family labor)
                 hfsup_bl(f,h)         factor supply by the household
                 fsup_bl(f)            factor supply
                 prev_bl(g,h)          production revenue from g by h
                 pcost_bl(g,h)         production cost for g by h
                 pprof_bl(g,h)         profit from activity g by h
                 finc_bl(f,h)          income from sale of factor f
                 efflabsup_bl          Wage-valued supply of effective labor at baseline
                 efflabsup_g_bl(g) Wage-valued effective labor used to produce g at baseline
                 labtime_g_bl(g)   Labor time supplied to produce g at baseline
                 totlabtime             Total amount of time supplied as labor (fixed over time)
                 harvest_ts(gstk,t)      Time series of harvest of stock good gstk
                 harvest_pc(gstk,t)      Percent change of harvest of stock good gstk

*infection dynamics parameters--baseline parameters (infectrate_BL is already defined above for model statements)
                beta_bl                  Baseline value of beta parameter
                gamma_bl                 Baseline value of gamma parameter
                lnagg_bl                 Baseline value of logged aggregate output (per capita)
                infectrate_BL      Baseline Infection Rate (endemic equilibrium)
                infectrate_hh_BL(h)      Baseline Infection Rate for household (endemic equilibrium)
                Ystate_BL                Baseline value of Y state variable
                fishlabtime_BL           Baseline value of total fishing labor time
                fishlabtime_hh_BL(h)     Baseline value of fishing labor time supplied by household
                vareps_hh_BL(h)          Baseline value of Household E-C Risk parameter (solved in baseline SS model)
                sigma_scalar_BL          Baseline value of R0-r0 wedge
                bigr0_bl                 Baseline value of R0
                flt_shr_BL(h)            Baseline value of Fishing labor time share for household
                extime_hh_bl(h)           Baseline value of total exposure time for household
;


*Define some fishing parameters:
Parameter
qp_fish                  value of quantity (millions) produced in fish sector--recall prices are 1 in baseline
harvest_kg1              Millions kg of fish harvested
aggr_price                     aggregate price of good (per kilo) for stocks (eg. fish) (derived from survey data)
kilo_conversion(g)             this factor converts value into kilos (kilos*UGX^-1)
X_frac_of_K                    assumption of what fraction standing stock biomass is of carrying capacity
d_sc(gstk,gg,h)           parameter for simple ratio search cost function
;

*===============================================================================
*================== DEFINE PARAMETERS FOR TIME-SERIES VALUES ===================
*=================== TO BE USED TO COMPUTE CHANGES OVER TIME ===================
*===============================================================================


* MACRO for defininf all parameters of type: ts, D, or PC

$macro defpars(i)  parameters    p&i(g,t)          price , \
         pva&i(g,h,t)       price of value added     ,   \
         qva&i(g,h,t)       quantity of value added  , \
         qp&i(g,h,t)            quantity produced by h        , \
         tqp&i(g,t)             total quantity produced of g  , \
         ttqp&i(t)             total production value in whole economy  , \
         ttqp_percap&i(t)       per capita production value in whole economy  , \
         Qp_composite&i(g,t)    Quntity of the composite good produced by armington func , \
         P_composite&i(g,t)     price of composite good                                 , \
         Imports&i(g,t)         imports used in the creation of the composite good      , \
         hqp&i(h,t)             total qty produced by a household  , \
         fd&i(g,f,h,t)          factor demand                           , \
         id&i(g,gg,h,t)         intermediate demand                      , \
         pshift&i(g,h,t)        cobb-douglas shifter                    , \
         fshare&i(g,f,h,t)      cobb-douglas shares                     , \
         r&i(g,f,h,t)           rent for fixed factors                  , \
         wv&i(f,t)            village-wide wage for tradable factors  , \
         vash&i(g,h,t)          value-added share                       , \
         idsh&i(gg,g,h,t)       intermediate demand share               , \
         tidsh&i(gg,h,t)        total intermediate input share (1-vash) , \
         fixfac&i(g,f,h,t)      fixed factor demand                     , \
         exinc&i(h,t)           exogenous income                        , \
         endow&i(f,h,t)         endowment                               , \
         qc&i(g,h,t)            level of consumption                    , \
         eshare&i(g,h,t)        consumption shares                      , \
         y&i(h,t)               income of household                     , \
         mry&i(t)                mean income in economy (weighted by nhh), \
         rytheil&i(t)            theil index for income                  , \
         cpi&i(h,t)             cpi                                     , \
         vqc&i(g,t)           village consumption                     , \
         vcpi&i(t)            village cpi                             , \
         cri&i(f,t)           rent weighted index                     , \
         ry&i(h,t)              real income                             , \
         ty&i(t)                income total                            , \
         try&i(t)               real income total                       , \
         emin&i(g,h,t)          incompressible demand                   , \
         trin&i(h,t)            transfers in - received                 , \
         trout&i(h,t)           transfers out - given                   , \
         sav&i(h,t)             savings                                 , \
         exproc&i(h,t)          expenditure rest of country             , \
         trinsh&i(h,t)          share of all transfers in the eco going to h  , \
         troutsh&i(h,t)         share of household h income being given as transfers , \
         hfd&i(f,h,t)           factor demand of household h for factor f , \
         vfd&i(f,t)           village demand for factor f             , \
         hms&i(g,h,t)           household marketed surplus of good g    , \
         vms&i(g,t)           village marketed surplus of good g       , \
         hfms&i(f,h,t)          household factor marketed surplus       , \
         vfms&i(f,t)          village factor marketed surplus         , \
         hfsup&i(f,h,t)         factor supply by the household                     , \
         fsup&i(f,t)            factor supply                     ,  \
         prev&i(g,h,t)          production revenue from g by h  , \
         pcost&i(g,h,t)         production cost for g by h      , \
         pprof&i(g,h,t)         profit from activity g by h    , \
         finc&i(f,h,t)          income from sale of factor f , \
         beta&i(t)               beta param    , \
         gamma&i(t)              gamma param   , \
         lnagg&i(t)              logged aggregate output (per capita) , \
         ystate&i(t)             Y state var                   , \
         infectrate&i(t)         Infection rate , \
         infectrate_hh&i(h,t)   Infection rate for household , \
         fishlabtime&i(t)        Total fishing labor time supplied by household , \
         fishlabtime_hh&i(h,t)  Total fishing labor time supplied by household,  \
         extime_ratio_hh&i(h,t)  Ratio of total exposure time for household at time t to baseline , \
         bigr0&i(t)              Value of R0 , \
         efflabsup&i(t)          Wage-valued supply of effective labor, \
         efflabsup_g&i(g,t)        Wage-valued effective labor used to produce g, \
         labtime_g&i(g,t)  Labor time supplied to produce g, \
         x&i(g,t)                stock level at time t (for use in fish diff eq) note this is indexed by good and time \
;

defpars(_ts) ;
defpars(_d) ;
defpars(_pc) ;


parameter theta_fish1_dr(g,f,h)       coefficient to represent share of the stocks value added
          sum_non_stock_betas(gstk,h) sum of the betas on the variable factors;

parameter y_dr3(h) income computed from total expenditures ;

parameter qcshare(h,g) share of household h in total consumption of g ;
parameter netexpsh(g) net export share of a good out of the zoi;
parameter qpshare(h,g) share of household h in production of g ;

parameter rough_id(g)                    temporary total intermediate demand
          share_dom_sector(g,gg)         the share of good gg that is dom prod and that is purch as ID in prod of g
          rough_tqp(g)                   rough estimate of total quantity produced for good g (whole economy)
          total_demanded(g)              total demand of good g found with sum of direct consumption and IDs
          total_dom_demanded(g)          total amount demanded that comes from domestic production sources
          share_dom(g)                   share of domestically produced goods in total amount of goods demanded
;

parameter fd_dr_check(g,f,h) Check calculation for fd_dr ;
parameter sumbetas_dr(g,h) Sum of production betas (inclusive of stock beta) ;

parameter qp_check(g,h) check QP with IRS fishing
          qva_check(g,h) check QVA with IRS fishing
;

parameter tid_dr(g) check of total id
          tqcid_dr(g)  check of qc+id
;

parameter
        shlab(h) share of village family labor coming from a household
;

parameter delta(g)                 share parameter for the armington CES aggregation function
          armg_shifter(g)          shift coefficient for the armington CES aggregation function
          p_import_dr(g)           the price of imports of good type g
          p_composite_dr(g)        price of the composite
          tqp_composite(g)         total qp of composite
          test_composite_output(g) test of composite output value;

parameter exinc_dr2(h) exogenous income computation
          exincsh2(h)  share of income being exogenous using exinc2
          feinc_dr(h)  income from factor endowments in the household
          fecomp_dr(f,h) income components ;

parameter        utility_sh(g,h)        CES share parameters for HHs
                 utility_coeff(h)       coefficient on CES utility function
                 qc_ratios(g,gg,h)      ratios of initial consumption of various goods for household h--g is numerator
                 sum_utility_sh(h)      test to make sure household share parameters add up to unity
                 calibr_demand(g,h)     test of demand from calibrated CES utility function--for verification
                 demand_test(g,h)       a check to see if calibrated CES demand is the same as the intial values of demand
                 p_dr(g)                prices (only for utility calibration below)
;


* =============================================================================
* ====================== CALIBRATE THE EPI-LEWIE ==============================
* ====================== COMPONENTS OF THE MODEL ==============================
* =============================================================================

*===================================================
*=========== LOCAL ECONOMY PARAMETERS ==============
*===================================================

* default values for parameters that were drawn
fshare_dr(g,f,h) = xlfshare(g,f,h) ;
fshare_dr(g,f,h)$fshare_dr(g,f,h) = fshare_dr(g,f,h)/sum(fa,fshare_dr(g,fa,h)) ;
eshare_dr(g,h)   = xleshare(g,h)   ;
troutsh_dr(h)   = xltroutsh(h)  ;
savsh_dr(h)     = xlSAVsh(h)    ;
exprocsh_dr(h)  = xlexpoutsh(h) ;

* quantity produced
qp_dr(g,h) = xlqp(g,h)*xlnhh(h) ;

*====================================================
*=========== FISH POP MODEL PARAMETERS ==============
*====================================================

*Add up the betas on non-stock factors:
sum_non_stock_betas(gstk,h) = sum(f,fshare_dr(gstk,f,h));
*Define and calculate each non-stock factors' share in the total value added
*across all non-stock factors:
theta_fish1_dr(g,f,h)$gp(g) = 0;
theta_fish1_dr(gstk,f,h)$fshare_dr(gstk,f,h) = fshare_dr(gstk,f,h) / sum_non_stock_betas(gstk,h);
display  theta_fish1_dr, sum_non_stock_betas, fshare_dr ;


*====================================================
*============= RETURN TO LE PARAMETERS ==============
*====================================================

p_dr(g) = 1 ;
r_dr(g,fk,h)     = 1 ;
wv_dr(ft)     = 1 ;


* Start from income

y_dr(h) = xlhhexp(h)*xlnhh(h) ;
display y_dr  ;

* all prices are 1 so cpi is 1
cpi_dr(h) = 1 ;
ry_dr(h) = y_dr(h) ;

* levels of expenditures on everything outside of the economy:
trout_dr(h) = y_dr(h)*troutsh_dr(h) ;
sav_dr(h) = y_dr(h)*savsh_dr(h) ;
exproc_dr(h) = y_dr(h)*exprocsh_dr(h) ;
display troutsh_dr, savsh_dr, exprocsh_dr, trout_dr, sav_dr, exproc_dr ;


* LEVELS OF CONSUMPTION:
display  p_dr,y_dr,sav_dr,trout_dr,exproc_dr, eshare_dr,g,h ;

qc_dr(g,h) = (y_dr(h)-sav_dr(h)-trout_dr(h)-exproc_dr(h))*eshare_dr(g,h)/p_dr(g) ;

display qc_dr ;

qcshare(h,g)$qc_dr(g,h) = qc_dr(g,h) / sum(hh,qc_dr(g,hh)) ;
display qcshare;

* PRODUCTION:

netexpsh(g)$gnag(g) = 1-(1/card(h)*(sum(h,xlrevsh_vil(g,h)+xlrevsh_zoi(g,h)))) ;
display netexpsh ;

* intermediate demand requirements
display xlQP;
idsh_dr(gg,g,h) = xlidsh(gg,g,h) ;
tidsh_dr(g,h) = sum(gg,idsh_dr(g,gg,h));
display idsh_dr, tidsh_dr;

tqc_dr(g) = sum(h,qc_dr(g,h)) ;
display tqc_dr ;

tqp_dr(g) = sum(h, qp_dr(g,h))  ;
ttqp_dr= sum(g,tqp_dr(g));
display tqp_dr, ttqp_dr ;

qpshare(h,g)$gp(g) = qp_dr(g,h)/sum(hh,qp_dr(g,hh)) ;
display qp_dr, qpshare ;


* all factor demands and intermediate demands:
id_dr(gg,g,h) = qp_dr(gg,h) * idsh_dr(gg,g,h) ;
display id_dr ;
display idsh_dr ;

*Create rough estimates of the share of each good's comsumption that comes from
*domestic production.
share_dom_sector(g,gg) = alldata("share_dom_sec",g,gg,"","FishPoor");

share_dom_sector(g,gg)$(share_dom_sector(g,gg)<.001) = 1 ;

rough_id(g) = sum((gg,h),qc_dr(gg,h)*(idsh_dr(gg,g,h)/(1-idsh_dr(gg,g,h)))) ;


*Rough estimate of total quantity produced
rough_tqp(g) = sum(h, qc_dr(g,h)) + rough_id(g) ;
*Total demanded
total_demanded(g) =  sum(h, qc_dr(g,h)) + sum(gg,rough_tqp(gg)*idsh_dr(gg,g,"FishPoor"));
*Total demanded from domestic sources
total_dom_demanded(g)  = sum(h, qc_dr(g,h)) + sum(gg,rough_tqp(gg)*idsh_dr(gg,g,"FishPoor")*share_dom_sector(gg,g));
*Calculate the share demanded that is domestically produced (condition on quantity consumed to exclude oil palm):
share_dom(g)$(tqc_dr(g)>0) = total_dom_demanded(g)/ total_demanded(g);

display rough_tqp, qc_dr, idsh_dr, share_dom_sector, total_demanded, total_dom_demanded, share_dom;

imports_dr(g)$(tqc_dr(g)>0) = ((1- share_dom(g))/share_dom(g))*tqp_dr(g);

* Factor demands derived from factor shares
fd_dr(g,f,h)$gp(g) = (qp_dr(g,h) - sum(gg,id_dr(g,gg,h))) * [fshare_dr(g,f,h)]  ;

fd_dr_check(g,f,h)$gp(g) = (qp_dr(g,h) - sum(gg,id_dr(g,gg,h))) * [fshare_dr(g,f,h) + theta_fish1_dr(g,f,h)*stockbeta(g)]  ;

fd_dr(g,f,h)$gp(g) = fd_dr_check(g,f,h) ;

*scale down factor demands in the case of IRS fishing production
sumbetas_dr(g,h)$gp(g) = sum(f,fshare_dr(g,f,h))+stockbeta(g)$(stockbeta(g) ne 1) ;
display sumbetas_dr ;
fd_dr(g,f,h)$gp(g) = fd_dr(g,f,h)/sumbetas_dr(g,h) ;

qva_dr(g,h)   = sum(f, fd_dr(g,f,h)) ;

pshift_dr(g,h)$(qva_dr(g,h))    = qva_dr(g,h)/prod(f,fd_dr(g,f,h)**fshare_dr(g,f,h)) ;

*compute value added share for all activities
vash_dr(g,h)$qp_dr(g,h) = (qp_dr(g,h)-sum(gg, id_dr(g,gg,h))) / qp_dr(g,h) ;
display id_dr, idsh_dr, tidsh_dr, vash_dr ;

qp_check(g,h)$qp_dr(g,h) = (1/vash_dr(g,h))*pshift_dr(g,h)*prod(f,fd_dr(g,f,h)**fshare_dr(g,f,h)) ;

qva_check(g,h)$qp_dr(g,h) = qp_check(g,h)*vash_dr(g,h) ;

display qp_dr, qp_check, qva_dr, qva_check ;

tid_dr(g)= sum((gg,h),id_dr(gg,g,h)) ;
tqcid_dr(g) = tid_dr(g) + tqc_dr(g) ;
display tqc_dr, tid_dr, tqcid_dr, tqp_dr ;


* FACTOR ENDOWMENTS :
* --------------------------
* for fixed factors, endowment is just factor use:
endow_dr(fk,h) = sum(g,fd_dr(g,fk,h)) ;
fixfac_dr(g,fk,h) = fd_dr(g,fk,h) ;

display endow_dr, fixfac_dr ;


* For labor
shlab(h) = (xlwrkagepop(h)*xlnhh(h))/sum((hh),xlwrkagepop(hh)*xlnhh(hh)) ;
shhh(h) = (xlnhh(h))/sum(hh,xlnhh(hh)) ;

display shlab, shhh ;

endow_dr("Labor",h) = shhh(h) * sum((hh,g), fd_dr(g,"Labor",hh)) ;
display endow_dr ;

display endow_dr ;


*Deriving parameters related to the Argminton function
*Normalize prices of imports and composite to 1:
p_import_dr(g) = 1;
p_composite_dr(g) = 1;

*Define the CES coefficient:
rho(g) = (trade_elas(g) - 1)/trade_elas(g);

*Calibrate the share parameter
delta(gm) = (imports_dr(gm)/tqp_dr(gm))**(rho(gm)-1) / ( (p_import_dr(gm)/p_dr(gm)) + (imports_dr(gm)/tqp_dr(gm))**(rho(gm)-1) );

*Set amount of composite good equal to the sum of domestic good & imports:
tqp_composite(g) = imports_dr(g) + tqp_dr(g);

*Calibrate the Armington shifter:
armg_shifter(gm) = tqp_composite(gm) /( delta(gm)*tqp_dr(gm)**(rho(gm)) + (1-delta(gm))*imports_dr(gm)**(rho(gm)))**(1/rho(gm)) ;

*Test of the output value created by the armington function to verify it
*is working.
test_composite_output(gm) = armg_shifter(gm) * ( delta(gm)*tqp_dr(gm)**(rho(gm)) + (1-delta(gm))*imports_dr(gm)**(rho(gm)))**(1/rho(gm));

display  delta,rho , armg_shifter,tqp_dr, imports_dr,tqp_composite, test_composite_output;



* MARKETS AGGREGATES
* ================================================================================================
* factor demand aggregates
hfd_dr(f,h)= sum(g,fd_dr(g,f,h)) ;
vfd_dr(f)= sum(h, hfd_dr(f,h)) ;

* marketed surpluses for goods
hms_dr(g,h) = qp_dr(g,h)$(not gm(g))  - qc_dr(g,h) - sum(gg,id_dr(gg,g,h)) ;
vms_dr(g) = tqp_composite(g)$gm(g) + sum(h,hms_dr(g,h)) ;

* marketed surpluses for factors
hfms_dr(ft,h) = endow_dr(ft,h) - sum(g, fd_dr(g,ft,h));
vfms_dr(ft) = sum(h, hfms_dr(ft,h));

* fixed factor demands at village
vfmsfix_dr(ftv) = vfms_dr(ftv) ;

* fixed goods trade levels at village
vmsfix_dr(gtv) = vms_dr(gtv) ;

* minimum consumption: zero for now.
emin_dr(g,h) = 0 ;

pva_dr(g,h) = p_dr(g) - sum(gg,idsh_dr(g,gg,h)*p_dr(gg)) ;

trinsh_dr(h)$(sum(hh,y_dr(hh)*xltrinsh(hh))) = y_dr(h)*xltrinsh(h)/sum(hh,y_dr(hh)*xltrinsh(hh))  ;
trin_dr(h) = trinsh_dr(h)*sum(hh,trout_dr(hh)) ;


* this is if we make exogenous income the residual from Y-FD
feinc_dr(h) = sum((g,fk),r_dr(g,fk,h)*fd_dr(g,fk,h)) + sum(ft, wv_dr(ft)*endow_dr(ft,h)) ;

* Make exogenous income / exogenous expenditures depending on what the sign is:
exinc_dr(h) = y_dr(h) - feinc_dr(h) ;


display pshift_dr, fshare_dr, p_dr, pva_dr, qva_dr, fd_dr, id_dr, r_dr, wv_dr, qp_dr, fixfac_dr, pva_dr,
        exinc_dr, endow_dr, y_dr,  qc_dr, eshare_dr, troutsh_dr, hfd_dr, vfd_dr,
        hms_dr, vms_dr, hfms_dr, vfms_dr ;


***Calibrating share parameters for CES utility function for households:
p_dr(g)=1;

*Define the  coefficient on CES utility function (not currently using this):
utility_coeff(h) = (good_elas(h) - 1) / good_elas(h);

*Define the consumption ratios
qc_ratios(g,gg,h)$(qc_dr(gg,h) > 0)=  qc_dr(g,h) / qc_dr(gg,h);

*Define the CES utility function share parameters:
utility_sh(g,h)$(qc_dr(g,h) > 0) = 1 / [ sum(gg, [qc_ratios(g,gg,h)**( -1/(good_elas(h)-1) )]$(qc_ratios(g,gg,h)>0) )];
*sum the shares up as a test that they add up to unity:
sum_utility_sh(h) = sum(g,utility_sh(g,h) );

*Now calculate demand using CES utility function and these utility_shares
calibr_demand(g,h) = [ p_dr(g)**(1-good_elas(h)) *  (utility_sh(g,h))**(good_elas(h)-1) *(y_dr(h)-sav_dr(h)-trout_dr(h)-exproc_dr(h))]$(utility_sh(g,h)>0) /
                     [ p_dr(g) * sum(gg, (p_dr(gg)**(1-good_elas(h)) *  utility_sh(gg,h)**(good_elas(h)-1))$(utility_sh(gg,h)>0) ) ] ;

*check to see if CES demand is the same as the intial values of demand:
*(These should be very close to zero, relative to scale of demand)
demand_test(g,h) =  calibr_demand(g,h) - qc_dr(g,h);

Display    qc_ratios, utility_sh, sum_utility_sh, calibr_demand, demand_test;


*=============================
*==== Epi Parameters =========
*=============================

alpha    = .15 ;
mu       = 0.017 ;

*beta
bmin     = 0.000580172388 ;
bmax     = 0.00580172388 ;
ymed     = 50 ;
yslope   = 10 ;

gmin     = 0.00055 ;
gmax     = 0.0055 ;

lnagg_dr = log(ttqp_dr*1000000/sum(h,xlnhh(h)*xlhhsize(h))) ;
*see stata do file "gdp distribution for gamma.do" for these next two values
lnagg_med  = 16.3791 ;
lnagg_s  = .54398003 ;

vareps_hh(h)       = alldata("Hhvareps","","","",h) ;
sigma_scalar        = sum(hh,shhh(hh)*vareps_hh(hh)**2) ;

display alpha, lnagg_dr, vareps_hh, sigma_scalar, lnagg_dr, lnagg_med, lnagg_s;


*==============================================================================
*==================== CALIBRATION OF BIOLOGICAL MODEL =========================
*==============================================================================

*Derive value of harvest
qp_fish = sum(h,qp_dr("fish",h));

*Aggregate fish price.
aggr_price("fish") = 5000;
*Make a coefficient to convert the value of harvest into kilos:
*Since qp_fish is in millions UGX, divide by the price per kilo to get millions of KG.
kilo_conversion("fish") = (1 /aggr_price("fish"));
*Derive harvest in kilograms:
harvest_kg1 =   qp_fish* kilo_conversion("fish") ;
*Save this as the baseline value of harvest:
harvestkg_bl("fish") =  harvest_kg1 ;

*Standing stock biomass is estimated to be 33% of carrying capacity(K):  X(t=0)=.33K
X_frac_of_K =  0.33 ;
*Now derive a carrying capacity such that growth equals kilos of harvest:
*If harvest = x(t)*gr*(1- x(t)/K), then
k_dr("fish") =  harvest_kg1 / [ (gr("fish")*X_frac_of_K*(1 - X_frac_of_K)) ];
*Now set initial condition for stock
stock_dr("fish")  = X_frac_of_K * k_dr("fish");
x_ts("fish",t) = X_frac_of_K * k_dr("fish");

*pass the value of the stock for the draw to the model parameter
stock("fish") = stock_dr("fish") ;

*Set this value to be the baseline stock level (x1) for use in
*calculating percent change later:
x_bl("fish")  =    stock_dr("fish");

display pshift_dr ;

*Now adjust the Cobb-Douglas shifter value to account for stock:
pshift_dr("fish",h)  =   pshift_dr("fish",h) / (stock_dr("fish")**stockbeta("fish")) ;

*pass this new shift parameter value to the model parameter
pshift("fish",h) = pshift_dr("fish",h) ;

**Calibrate d_sc:
d_sc(gstk,gg,h) = (stock_dr(gstk)**n_sc)*idsh_dr(gstk,gg,h);

display  qp_fish, kilo_conversion, harvestkg_bl, harvest_kg1, gr, k_dr, stock_dr, x_ts, pshift_dr, aggr_price, n_sc, d_sc, idsh_dr  ;


*==============================================================================
*==================== SOLVING THE BASELINE MODEL  =============================
*==============================================================================


*===================================================================================
*=================== PREPARE TO SOLVE THE BASELINE MODEL  ==========================
*===================================================================================

* re-initialise all the variables in the matrix
pshift(g,h)    = pshift_dr(g,h) ;
fshare(g,f,h)  = fshare_dr(g,f,h) ;
P.l(g)        = p_dr(g) ;
QVA.l(g,h)     = qva_dr(g,h) ;
FD.l(g,f,h)    = fd_dr(g,f,h) ;
ID.l(gg,g,h)   = id_dr(gg,g,h) ;
R.l(g,fk,h)    = r_dr(g,fk,h);
WV.l(f)      = wv_dr(f) ;
QP.l(g,h)      = qp_dr(g,h) ;
fixfac(g,fk,h) = fixfac_dr(g,fk,h) ;
vfmsfix(ftv) = vfmsfix_dr(ftv) ;
PVA.l(g,h)     = pva_dr(g,h) ;
vash(g,h)      = vash_dr(g,h) ;
idsh(gg,g,h)   = idsh_dr(gg,g,h) ;
tidsh(g,h)     = tidsh_dr(g,h) ;
exinc(h)       = exinc_dr(h) ;
endow(f,h)     = endow_dr(f,h) ;
Y.l(h)         = y_dr(h) ;
CPI.l(h)       = cpi_dr(h) ;
RY.l(h)        = ry_dr(h) ;
TRIN.l(h)      = trin_dr(h) ;
trinsh(h)      = trinsh_dr(h) ;
emin(g,h)      = emin_dr(g,h) ;
QC.l(g,h)      = qc_dr(g,h) ;
eshare(g,h)    = eshare_dr(g,h) ;
troutsh(h)     = troutsh_dr(h) ;
TROUT.l(h)     = trout_dr(h) ;
HFD.l(f,h)     = hfd_dr(f,h);
VFD.l(f)     = vfd_dr(f);
HMS.l(g,h)     = hms_dr(g,h);
VMS.l(g)     = vms_dr(g);
vmsfix(gtv)  = vmsfix_dr(gtv);
HFMS.l(ft,h)   = hfms_dr(ft,h);
VFMS.l(ft)   = vfms_dr(ft);
savsh(h)       = savsh_dr(h) ;
exprocsh(h)    = exprocsh_dr(h) ;
SAV.l(h)       = sav_dr(h) ;
EXPROC.l(h)    = exproc_dr(h) ;
hfsupzero(ft,h) = endow_dr(ft,h) ;

utility_sh_m(g,h) = utility_sh(g,h) ;
delta_m(g)        = delta(g) ;
armg_shifter_m(g) = armg_shifter(g) ;
p_import_m(g)     = p_import_dr(g) ;

lnagg_m           = lnagg_dr  ;

sumbetas(g,h)  = sumbetas_dr(g,h) ;

*Variables: related to the Armington function:
QP_COMPOSITE.l(gm) = tqp_composite(gm) ;
P_COMPOSITE.l(gm)  = p_composite_dr(gm);
IMPORTS.l(gm)      = imports_dr(gm)  ;

* Make sure FD cannot reach zero if it wasn't at zero:
FD.lo(g,f,h)$FD.l(g,f,h) = 1E-15;

*fish stock
stock(g) = stock_dr(g) ;
theta_fish1(g,f,h) = theta_fish1_dr(g,f,h) ;
sumbetas(g,h)    = sumbetas_dr(g,h) ;

*set initial guesses for Epi variables
fishlabtime_hh.l(h) = fd_dr("fish","labor",h)/(1-0.5*alpha) ;
fishlabtime.l      = sum(h,fishlabtime_hh.l(h)) ;
beta.l             = (1 + exp((lnagg_dr-lnagg_med)/lnagg_s))**(-1)*(bmax-bmin) + bmin ;
gamma.l             = (1 + exp(-(lnagg_dr-lnagg_med)/lnagg_s))**(-1)*(gmax-gmin) + gmin ;
Infectrate_hh.l(h)  = 0.5 ;
Ystate.l           = 0.5 ;

*** Labor Supply variables
hfsupel("LABOR",h) = 100 ;

*hfsupel("FL",h) = %flse% ;
HFSUP.l(f,h)    = hfsupzero(f,h) ;
display HFSUP.l, hfsupzero ;

*** Baseline labor time
LABTIME_BLvar.l(f,h) =  sum(g,fd_dr(g,f,h))/(1-0.5*alpha) ;


* closures: fixed wages and prices on world-market-integrated factors and goods (ftw & gtw)
WV.fx(ftw) = WV.l(ftw);
P.fx(gtw) = P.l(gtw) ;

*this is to fix the problem that this was an unmatched variable in the baseline solve.
*nobody produces OUT, nobody consumes palmoil, and nonfishing households dont produce fish
FD.fx("OUT",f,h)=0 ;
QP.fx("OUT",h)=0 ;
QC.fx("palmoil",h)=0 ;
QP.fx("fish","NfishNpoor")=0 ;
QP.fx("fish","NfishPoor")=0 ;

*---------------------------------
* SOLVING THE BASELINE MODEL
*---------------------------------
option iterlim = 100000 ;

BL_MODEL.limcol = 0 ;
BL_MODEL.limrow = 0 ;
solve BL_MODEL using mcp ;
ABORT$(BL_MODEL.modelstat ne 1) "NOT WELL CALIBRATED IN THIS DRAW - CHECK THE DATA INPUTS" ;
display BL_MODEL.modelstat ;

pshift_bl(g,h)   = pshift(g,h) ;
fshare_bl(g,f,h) = fshare(g,f,h) ;
p_bl(g)       = P.l(g) ;
qva_bl(g,h)      = QVA.l(g,h) ;
fd_bl(g,f,h)     = FD.l(g,f,h) ;
id_bl(gg,g,h)    = ID.l(gg,g,h) ;
r_bl(g,fk,h)     = R.l(g,fk,h) ;
wv_bl(f)       = WV.l(f) ;
vash_bl(g,h)     = vash(g,h)  ;
qp_bl(g,h)       = QP.l(g,h) ;
fixfac_bl(g,fk,h) = fixfac(g,fk,h) ;
pva_bl(g,h)      = PVA.l(g,h) ;
vash_bl(g,h)     = vash(g,h) ;
idsh_bl(g,gg,h)  = idsh(g,gg,h) ;
tidsh_bl(g,h)    = tidsh(g,h) ;
exinc_bl(h)      = exinc(h) ;
endow_bl(f,h)    = endow(f,h) ;
y_bl(h)          = Y.l(h) ;
qc_bl(g,h)       = QC.l(g,h) ;
cpi_bl(h)        = CPI.l(h) ;
vqc_bl(g)      = sum(h, qc_bl(g,h));
* village cpi is weighted sum of prices
vcpi_bl       = sum((h,g), (p_bl(g)**2)*qc_bl(g,h)) / sum((h,g),p_bl(g)*qc_bl(g,h)) ;
cri_bl(f)      = sum((g,h), r_bl(g,f,h)*fd_bl(g,f,h)/sum((gg,hh),fd_bl(gg,f,hh)) ) ;
ry_bl(h)         = RY.l(h) ;
ty_bl            = sum(h,y_bl(h));
try_bl           = sum(h,ry_bl(h));
trin_bl(h)       = TRIN.l(h) ;
trout_bl(h)      = TROUT.l(h) ;
trinsh_bl(h)     = trinsh(h) ;
sav_bl(h)        = SAV.l(h) ;
exproc_bl(h)     = EXPROC.l(h) ;
eshare_bl(g,h)   = eshare(g,h) ;
emin_bl(g,h)     = emin(g,h) ;
troutsh_bl(h)    = troutsh(h) ;
hfd_bl(f,h)      = HFD.l(f,h) ;
vfd_bl(f)      = VFD.l(f) ;
hms_bl(g,h)      = HMS.l(g,h) ;
vms_bl(g)      = VMS.l(g) ;
hfms_bl(ft,h)    = HFMS.l(ft,h) ;
vfms_bl(ft)    = VFMS.l(ft) ;
hfsup_bl(ft,h)   = HFSUP.l(ft,h) ;
fsup_bl(ft)      = sum(h,hfsup_bl(ft,h)) ;

Qp_composite_bl(g) = Qp_composite.l(g) ;
P_composite_bl(g)  = P_composite.l(g) ;
Imports_bl(g)      = Imports.l(g) ;

* Factor income from sales
finc_bl(f,h)      = HFMS.l(f,h)*(WV.l(f)$(ftw(f)) + WV.l(f)$ftv(f)) ;

vfmsfix_bl(ft) = vfmsfix_dr(ft) ;

* more params
tqp_bl(g)        = sum(h,qp_bl(g,h)) ;
ttqp_bl         = sum(g,tqp_bl(g)) ;
ttqp_percap_bl  = ttqp_bl / sum(h,xlnhh(h)*xlhhsize(h)) ;
hqp_bl(h)        = sum(g, qp_bl(g,h)) ;
efflabsup_bl     = sum(h, sum(g, fd_dr(g,"labor",h))) ;
efflabsup_g_bl(g) = sum(h, fd_dr(g,"labor",h)) ;
labtime_g_bl(g)  = sum(h, fd_dr(g,"labor",h)/(1-Infectrate_hh.l(h)*alpha)) ;

* mean income and theil index:
mry_bl           = sum(h, ry_bl(h)) / sum(h,xlnhh(h)) ;
* modified formula a bit to account for the fact that we grouped the households
rytheil_bl        = sum(h, (ry_bl(h)/xlnhh(h))/mry_bl
                             *log((ry_bl(h)/xlnhh(h))/mry_bl)
                             *xlnhh(h))
                          / sum(h,xlnhh(h)) ;

* Pass the solution for baseline labor time to the parameter for the SIMS model
labtime_bl(f,h)          = LABTIME_BLvar.l(f,h) ;
totlabtime       = sum(h,LABTIME_BLvar.l("labor",h)) ;


*Pass the model solutions to the baseline parameters
fishlabtime_hh_BL(h)     = fishlabtime_hh.l(h) ;
fishlabtime_BL           = fishlabtime.l ;


extime_hh_bl(h)             = (fishlabtime_hh_BL(h)*1.25)$fishlabtime_hh_BL(h) + 1$(not fishlabtime_hh_BL(h)) ;
*extime_hh_bl(h)             = 1 ;

Infectrate_hh_BL(h)             = Infectrate_hh.l(h) ;
infectrate_BL = sum(h,infectrate_hh_BL(h)*shhh(h)) ;

Ystate_BL                = Ystate.l ;
beta_BL                  = beta.l ;
gamma_BL                 = gamma.l ;
lnagg_bl                 = (sum(g,sum(h,QP.l(g,h)))*1000000/sum(h,xlnhh(h)*xlhhsize(h))) ;


* Pass baseline values to 1st period values for infectrate_hh_ts and Ystate_ts
Ystate_ts("1") = Ystate_BL ;
infectrate_hh_ts(h,"1") = infectrate_hh_BL(h) ;

*Pass 1st period values for infectrate_hh to sims parameter in SIMS model
Infectrate_hh_sims(h) = infectrate_hh_ts(h,"1") ;

*calculate R0
bigr0_bl = beta_BL*beta_BL*sum(h,shhh(h)*vareps_hh(h)**2)/(gamma_bl*mu) ;

display stock, stock_dr, pshift, pshift_dr, ttqp_bl, tqp_bl;


*======================================================
*================ ONE-TIME SHOCKS =====================
*======================================================


* Treatment that increases mortality of parasite in human host
*===============================================================


*======================================================
* 25% restriction on fishing capital
*======================================================
fixfac("fish","capital",h)$(FMPshock = 1) = fixfac("fish","capital",h)*(1-25/100) ;
FD.l("fish","capital",h)$(FMPshock = 1) = fixfac("fish","capital",h) ;


*=============================================================================
*=============== TIME SERIES LOOP BEGINS HERE ================================
*=============================================================================

*set counter =1 for if/else statement in fish transition equation
tt=1 ;

loop (t,

*Adjust the idsh for the fishing sector
idsh(gstk,gg,h)= d_sc(gstk,gg,h)/ (stock(gstk)**n_sc);

*====================================================
*================ ANNUAL SHOCKS =====================
*====================================================

*1% increase in Oil Palm TFP each year
*=====================================
pshift("palmoil",h)$(TFPshock = 1) = pshift("palmoil",h)*1.01 ;

* Using PATH solver and MCP model
SIMS_MODEL.limcol = 0 ;
SIMS_MODEL.limrow = 0 ;

SIMS_MODEL.solveOpt=0 ;

solve SIMS_MODEL using mcp ;
ABORT$(SIMS_MODEL.modelstat ne 1) "NO OPTIMAL SOLUTION REACHED in MCP model" ;
solve SIMS_MODEL using mcp ;
modstat = SIMS_MODEL.modelstat ;
modstat_dr = SIMS_MODEL.modelstat ;

*==================================================
*====== UPDATE BIOLOGICAL MODEL PARAMETERS ========
*==================================================

*Use the production in the fishing sector to calculate harvest in kilograms:
         harvestkg_ts(gstk,t) =  sum(h,QP.l(gstk,h))*kilo_conversion(gstk) ;
*Derive fish population size for end of time period t

if ((tt<2),
    x_ts(gstk,t) = x_bl(gstk) + gr(gstk)*x_bl(gstk)*(1-x_bl(gstk)/k_dr(gstk)) - sum(h,QP.l(gstk,h))*kilo_conversion(gstk) ;
else
    x_ts(gstk,t) = x_ts(gstk,t-1) + gr(gstk)*x_ts(gstk,t-1)*(1-x_ts(gstk,t-1)/k_dr(gstk)) - sum(h,QP.l(gstk,h))*kilo_conversion(gstk) ;
);

*Redefine the model stock 'parameter' for use in t+1 LEWIE.
stock(gstk) = x_ts(gstk,t);

*===================================================================
*====== SAVE OUTPUT FROM EACH LOOP IN TIME-SERIES VARIABLES ========
*===================================================================

pshift_ts(g,h,t)   = pshift(g,h) ;
fshare_ts(g,f,h,t) = fshare(g,f,h) ;
p_ts(g,t)       = P.l(g) ;
qva_ts(g,h,t)      = QVA.l(g,h) ;
fd_ts(g,f,h,t)     = FD.l(g,f,h) ;
id_ts(gg,g,h,t)    = ID.l(gg,g,h) ;
r_ts(g,fk,h,t)     = R.l(g,fk,h) ;
wv_ts(f,t)       = WV.l(f) ;
idsh_ts(gg,g,h,t)  = idsh(gg,g,h) ;
tidsh_ts(g,h,t)    = tidsh(g,h) ;
qp_ts(g,h,t)       = QP.l(g,h) ;
tqp_ts(g,t)        = sum(h,qp_ts(g,h,t)) ;
ttqp_ts(t)         = sum(g,tqp_ts(g,t)) ;
ttqp_percap_ts(t)  = ttqp_ts(t) / sum(h,xlnhh(h)*xlhhsize(h)) ;
hqp_ts(h,t)        = sum(g, qp_ts(g,h,t)) ;

fixfac_ts(g,fk,h,t) = fixfac(g,fk,h) ;
pva_ts(g,h,t)      = PVA.l(g,h) ;
vash_ts(g,h,t)     = vash(g,h)  ;
exinc_ts(h,t)      = exinc(h) ;
endow_ts(f,h,t)    = endow(f,h) ;
y_ts(h,t)          = Y.l(h) ;
qc_ts(g,h,t)       = QC.l(g,h) ;
cpi_ts(h,t)        = CPI.l(h) ;
vqc_ts(g,t)      = sum(h, qc_ts(g,h,t));
* village cpi is weighted sum of prices
vcpi_ts(t)       = sum((h,g), (p_ts(g,t)**2)*qc_ts(g,h,t)) / sum((h,g),p_ts(g,t)*qc_ts(g,h,t)) ;
* weighted capital rent in the village
cri_ts(f,t)          = sum((g,h), r_ts(g,f,h,t)*fd_ts(g,f,h,t)/sum((gg,hh),fd_ts(gg,f,hh,t)) ) ;
ry_ts(h,t)         = RY.l(h) ;
ty_ts(t)           = sum(h,y_ts(h,t));
try_ts(t)          = sum(h,ry_ts(h,t));
trinsh_ts(h,t)     = trinsh(h) ;
eshare_ts(g,h,t)   = eshare(g,h) ;
emin_ts(g,h,t)     = emin(g,h) ;
troutsh_ts(h,t)    = troutsh(h) ;
hfd_ts(f,h,t)      = HFD.l(f,h) ;
vfd_ts(f,t)      = VFD.l(f) ;
hms_ts(g,h,t)      = HMS.l(g,h) ;
vms_ts(g,t)      = VMS.l(g) ;
hfms_ts(ft,h,t)    = HFMS.l(ft,h) ;
vfms_ts(ft,t)    = VFMS.l(ft) ;
trin_ts(h,t)       = TRIN.l(h) ;
trout_ts(h,t)      = TROUT.l(h) ;
sav_ts(h,t)        = SAV.l(h) ;
exproc_ts(h,t)     = EXPROC.l(h) ;
hfsup_ts(ft,h,t)   = HFSUP.l(ft,h) ;
fsup_ts(ft,t)      = sum(h,hfsup_ts(ft,h,t)) ;
efflabsup_ts(t)    = sum(h, sum(g, fd_ts(g,"labor",h,t))) ;
efflabsup_g_ts(g,t) = sum(h, fd_ts(g,"labor",h,t)) ;
labtime_g_ts(g,t)  = sum(h, fd_ts(g,"labor",h,t)/(1-Infectrate_hh_sims(h)*alpha)) ;

Qp_composite_ts(g,t) = Qp_composite.l(g) ;
P_composite_ts(g,t)  = P_composite.l(g) ;
Imports_ts(g,t)      = Imports.l(g) ;

idsh_fish_record_ts(gstk,gg,h,t)= idsh(gstk,gg,h) ;

* Factor income from sales
finc_ts(f,h,t)      = HFMS.l(f,h)*(WV.l(f)$(ftw(f)) + WV.l(f)$ftv(f)) ;


mry_ts(t)           = sum(h, ry_ts(h,t)) / sum(h,xlnhh(h)) ;
* modified formula a bit to account for the fact that we grouped the households
rytheil_ts(t)        = sum(h, ry_ts(h,t)/xlnhh(h)/mry_ts(t)
                             *log(ry_ts(h,t)/xlnhh(h)/mry_ts(t))
                             *xlnhh(h))
                          / sum(h,xlnhh(h)) ;

*===============================
*===== INFECTION DYNAMICS ======
*===============================

fishlabtime_hh_ts(h,t) = fishlabtime_hh.l(h) ;
fishlabtime_ts(t) = fishlabtime.l ;
lnagg_ts(t) = log(ttqp_ts(t)*1000000/sum(h,xlnhh(h)*xlhhsize(h))) ;
beta_ts(t) = beta.l ;
gamma_ts(t) = gamma.l ;
extime_ratio_hh_ts(h,t) = (0.2 + fishlabtime_hh.l(h)/extime_hh_bl(h))$fishlabtime_hh_ts(h,t) + 1$(not fishlabtime_hh_ts(h,t)) ;

*update mdaeffcov annually if running MDA shock. base of 0.295 provides a 19% reduction in I_h (King 2020)
mdaeffcov$(MDAshock = 1) = gamma_bl*0.295**(1/tt);

*Calculate R0
bigr0_ts(t)      = beta_ts(t)*beta_ts(t)*sum(h,shhh(h)*extime_ratio_hh_ts(h,t)**2*vareps_hh(h)**2)/((mdaeffcov+gamma_ts(t))*mu) ;


*Solve the Quasi-SS model
*First, *set initial guesses for variables using previous period values
         Infectrate_hh_ts_t.l(h)  = Infectrate_hh.l(h) ;
         Ystate_ts_t.l           = Ystate.l ;
*assign beta and gamma and extime ratio values calculated above to parameters in model
         beta_ts_t       =       beta_ts(t) ;
         gamma_ts_t      =       gamma_ts(t) ;
         extime_ratio_t(h)  =       extime_ratio_hh_ts(h,t) ;

         qssa.limcol = 0;
         qssa.limrow = 0;
         solve qssa using mcp ;
         ABORT$(qssa.modelstat ne 1) "NO OPTIMAL SOLUTION REACHED in MCP model" ;
         solve qssa using mcp ;

*Now pass on values from the solution to the state variables in time t. These are the EE values in time t
         infectrate_hh_ts(h,t) = Infectrate_hh_ts_t.l(h) ;
         ystate_ts(t) = Ystate_ts_t.l ;
*Pass t-period values for infectrate_hh to sims parameter in SIMS model
         Infectrate_hh_sims(h) = Infectrate_hh_ts_t.l(h) ;


*update counter
tt = tt+1 ;

***TIME LOOP ENDS HERE
display "THIS IS THE END OF THE TIME LOOP" ;
);

display fd_bl, fd_ts, fixfac_bl, fixfac_ts, bigr0_bl, bigr0_ts,
         infectrate_hh_bl, infectrate_hh_ts, infectrate_bl,
         fishlabtime_bl, fishlabtime_ts, fishlabtime_hh_bl, fishlabtime_hh_ts,
         ttqp_bl, tqp_ts, ttqp_ts, ry_bl, ry_ts,
         beta_bl, beta_ts, gamma_bl, gamma_ts, vareps_hh, ystate_bl, ystate_ts,
         sigma_scalar, efflabsup_bl, efflabsup_ts, p_ts, P_composite_ts, x_ts, idsh_fish_record_ts ;


execute_unload "EEL-output.gdx" ;
