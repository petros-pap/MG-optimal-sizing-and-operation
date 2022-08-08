sets t  time interval (one hour) /t1*t8760/

     mm maximum number of microturbines considered /m1, m2, m3/

     eq equipment considered in the model

  /  m  "Microturbines"
     w  "Wind turbines"
     s  "Solar PVs"
     j  "Electric boilers"
     n  "Natural gas boilers"
     b  "Battery bank"  /;

*~~~~~~~~~~~~~~~~~~~~~~~~~~~ Excel Data Extraction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set  sd statistical data of parameters /fw, H, Pl, Ql/;

*Where: fw = Availability fraction of Rated Power Output of wind turbines
*       H  = Hourly solar radiation
*       Pl = Power consumed by the Load
*       Ql = Heat consumed by the Load

* Conversion of Excel file to GDX file
$CALL GDXXRW.EXE GAMS_DATA_1.xlsx par=dataset rng=Sheet1!A1:E8761

*Iporting DATA from GDX to compilation script
parameter dataset(t,sd);
$GDXIN GAMS_DATA_1.gdx
$LOAD dataset
$GDXIN
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parameters

C(eq)            capital cost coefficient of the equipment

              /   m = 594000
                  w = 3400
                  s = 5000
                  j = 60
                  n = 60
                  b = 132      /

O(eq)            maintenance cost coefficient of the equipment

              /   m = 0.02
                  w = 0.008
                  s = 52
                  j = 0.0075
                  n = 0.0075
                  b = 0.00143  /;


scalars

Agrid            purchase cost of power from the macrogrid [$ per kWh]   /0.108/
Ang              purchase cost of natural gas [$ per kWh]                /0.0302/
gamma            heat-to-power ratio of microturbines [-]                /1.5/
Hr               reference solar radiation of PVs [kW per m2]            /1/
nm               electrical efficiency of microturbines [-]              /0.27/
nbm              charging efficiency of battery bank [-]                 /0.9/
nbp              discharging efficiency of battery bank [-]              /0.95/
nj               power conversion efficiency of the electric boiler [-]  /0.9/
nn               natural gas conversion efficiency of the NG boiler [-]  /0.85/
Ostartup         microturbine additional cost per startup [$ per start]  /10/
Pmr              rated power output for each microturbine [kW]           /165/
Hng              assumed heating value for natural gas [kWh per m3]      /9.83/
theta_n          NPV factor of capital cost of NG boiler                 /1.2792/
theta_j          NPV factor of capital cost of electric boiler           /1.2792/
theta_b          NPV factor of capital cost of battery bank              /2.4241/
phi_ng           NPV factor of annual purchase cost of natural gas       /11.7948/
phi_grid         NPV factor of annual purchase cost of power from grid   /11.7948/
phi_maint        NPV factor of maintenance cost of all the equipment     /11.7948/
Fng              natural gas CO2 emissions coefficient [tonCO2 per kWh]  /0.000183/
Fgrid            macrogrid CO2 emissions coefficient [tonCO2 per kWh]    /0.000575/;



binary variables

x(t,mm)          micriturbine mm is on
y(mm)            microturbine mm is installed
z(t,mm)          microturbine mm is started up;


positive variables

Capex            capital expenditures [USD]
Opex             operational expenditures [USD]
Pwr              rated power of windturbines [kW]
Psr              rated power of solar PVs [kW]
Qjr              rated heat of electric boiler [kW]
Qnr              rated heat of natural gas boiler [kW]
Lbmax            maximum energy storage of battery bank [kW]
Pm(t,mm)         microturbine hourly power output [kWh]
Pw(t)            windturbine power output [kWh]
Qn(t)            natural gas boiler hourly heat output [kWh]
Qj(t)            electric boiler hourly heat output [kWh]
Pbplus(t)        battery bank hourly discharge [kWh]
Gm(t,mm)         natural gas hourly consumption by the MTs [m3]
Gn(t)            natural gas hourly consumtion by the NGB [m3]
Pgrid(t)         hourly power purchased from the grid [kWh]
Ps(t)            hourly power produced by solar PVs [kWh]
Pbminus(t)       battery bank hourly charge [kWh]
Pj(t)            hourly electric boiler power consumption [kWh]
Pd(t)            power dumped each hour [kWh]
Qm(t,mm)         hourly heat produced by the MTs [kWh]
Qd(t)            dump heat load during an hour [kWh]
Lb(t)            hourly energy storage level [kWh]
En               gas fired boiler emissions [ton(CO2)]
Em               microturbine emissions [ton(CO2)]
Egrid            macrogrid emissions [ton(CO2)];

variable NPC     net present value of costs [USD];


equations

eq2              sum of the capital costs
eq3              sum of the annual operational costs
eq4(t)           electric power balance
eq5(t)           heat balance
eq11a(t,mm)      lower bound from MTs efficiency constraint
eq11b(t,mm)      upper bound from MTs efficiency constraint
eq12(t,mm)       useful heat generated from MTs
eq13(t,mm)       microtubines fuel consumption
eq14(t,mm)       number of activated MTs limitation by the number intalled
eq15(t,mm)       track of startups of microturbines
eq16(t)          solar PVs power balance
eq17(t)          wind turbine power balance
*eq20a(t)         lower bound of EB hourly consumption constraint
eq20b(t)         upper bound of EB hourly consumption constraint
*eq21a(t)         lower bound of NGB hourly consumption constraint
eq21b(t)         upper bound of NGB hourly consumption constraint
eq22(t)          heat to power conversion of EB
eq23(t)          natural gas conversion of NGB
eq25(t)          battery bank power balance
eq26a(t)         lower bound in battery bank storage level constraint
eq26b(t)         upper bound in battery bank storage level constraint
mtCO2            microturbine CO2 emissions
ngbCO2           gas fired boiler CO2 emissions
gridCO2          macrogrid CO2 emissions

obj              sum of capital and operational expenditures;


*~~~~~~~~~~~~~~~~ Annual Operational and Capital expenses ~~~~~~~~~~~~~~~~~~~~~~

eq2..            Capex =e= sum(mm,y(mm))*C('m')
                           +Pwr*C('w')
                           +Psr*C('s')
                           +Qjr*C('j')*theta_j
                           +Qnr*C('n')*theta_n
                           +Lbmax*C('b')*theta_b;

eq3..            Opex  =e= Psr*O('s')*phi_maint+sum(t,
                           sum(mm,(Pm(t,mm)*O('m')+z(t,mm)*Ostartup))*phi_maint
                           +Pw(t)*O('w')*phi_maint
                           +Qn(t)*O('n')*phi_maint
                           +Qj(t)*O('j')*phi_maint
                           +Pbplus(t)*O('b')*phi_maint
                           +(sum(mm,Gm(t,mm))+Gn(t))*Ang*phi_ng
                           +Pgrid(t)*Agrid*phi_grid );

*~~~~~~~~~~~~~~~~~~~~~~~~~~ Power and Heat Balance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eq4(t)..         sum(mm,Pm(t,mm))+Ps(t)+Pw(t)+Pbplus(t)+Pgrid(t) =e=
                 dataset(t,'Pl')+Pbminus(t)+Pj(t)+Pd(t);

eq5(t)..         sum(mm,Qm(t,mm))+Qj(t)+Qn(t) =e= dataset(t,'Ql')+Qd(t);

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Equipment Design ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Microturbines

eq11a(t,mm)..    0.5*x(t,mm)*Pmr =l= Pm(t,mm);

eq11b(t,mm)..    Pm(t,mm) =l= x(t,mm)*Pmr;

eq12(t,mm)..     Qm(t,mm) =e= Pm(t,mm)*gamma;

eq13(t,mm)..     Gm(t,mm) =e= Pm(t,mm)/nm;

eq14(t,mm)..     x(t,mm) =l= y(mm);

eq15(t,mm)$(ord(t) le card(t)-1)..
                 z(t+1,mm) =g= x(t+1,mm)-x(t,mm);

*Solar PV Array

eq16(t)..        Ps(t) =e= Psr*dataset(t,'H')/Hr;

*Wind Turbines

eq17(t)..        Pw(t) =e= dataset(t,'fw')*Pwr;

*Boilers

*eq20a(t)..       0 =l= Qj(t);

eq20b(t)..       Qj(t) =l= Qjr;

*eq21a(t)..       0 =l= Qn(t);

eq21b(t)..       Qn(t) =l= Qnr;

eq22(t)..        Pj(t) =e= Qj(t)/nj;

eq23(t)..        Gn(t) =e= Qn(t)/nn;

*Battery Bank

eq25(t)..        Lb(t) =e= Lb(t-1)+Pbminus(t)*nbm-Pbplus(t)/nbp;

eq26a(t)..       0.2*Lbmax =l= Lb(t);

eq26b(t)..       Lb(t) =l= Lbmax;

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CO2 Emisiions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mtCO2..          Em =e= sum(t,sum(mm,Gm(t,mm)))*Fng;

ngbCO2..         En =e= sum(t,Gn(t))*Fng;

gridCO2..        Egrid =e= sum(t,Pgrid(t))*Fgrid;

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Objective Function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

obj..            NPC =e= Opex + Capex;


model    Microgrid       /ALL/;

options  MIP   = CPLEX
         optcr = 0.01 ;

Microgrid.iterlim=1000000;

solve Microgrid using MIP minimizing NPC




*~~~~~~~~~~~~ Unload to GDX file (occurs during execution phase) ~~~~~~~~~~~~~~~

execute_unload "Ref_Case.gdx"    Pgrid.L Qn.L Qj.L Pw.L Ps.L Pm.L Qm.L Lb.L
                                 Capex.L Opex.L NPC.L y.L Psr.L Pwr.L Lbmax.L
                                 Qjr.L Qnr.L Em.L En.L Egrid.L;


$ONTEXT
file out /testgrid.dat/;
put out;

loop(t,
put "Time interval", ord(t):1:0, "Pgrid = " , Pgrid.L(t) ;
    );
putclose out;

Variable Extraction to Excel File

Open the GDX file created above and manually create an Excel file by
right clicking on any symbol and selecting 'Write all symbols to XLS file'
$OFFTEXT