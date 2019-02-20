clear all
close all
%run this script first
start_time = tic;

%Automation of runs setup:


%Properties for runs:
POWER_SA = 1400; %MW, maximum electric power provided by steam accumulator
T_STORE = 1*3600; %s, storage time
LossRateTarget = 0.1; %[%]
InsulationThickness = 1; %[m]


%SA properties
T_MAX = 24; %hr, time at maximum SA power
%LTANK = 500000*T_MAX/4 % m, pipe length
LTANK = 3150000;
D_RAMP_RATE = 1.67; %percent/min, discharge ramp rate (SA turbine, % of max power/min)
POWER_SA_INITIAL = 0; %MW, cold start
T_RAMP = 60*((POWER_SA-POWER_SA_INITIAL)/((D_RAMP_RATE/100)*POWER_SA)); %s, time to ramp up to max power
T_END = T_RAMP+T_MAX*3600; %s, discharge time including ramp up and time at max power
ENERGY_SA = 0.5*T_RAMP*(POWER_SA-POWER_SA_INITIAL)/3600+POWER_SA*T_MAX; %MWh, electric energy provided by SA in one discharge cycle
DT = 10;  % s, time step
RTANK = 0.4064; % m (16 inches)
VTANK = LTANK.*pi.*RTANK^2.;  % m3
MaxTemp = 285; %Used in loss rate/ hr calculation

t_length = floor(T_RAMP/DT); %length of ramping power vector
pow = zeros(t_length,1);
for t = 1:t_length
    pow(t) = t*DT*((D_RAMP_RATE/100)*POWER_SA)/60; %MW, power as SA turbine is ramping up
end
%accumulator initial thermo properties
P0 = 70; % bar, initial pressure
X0 = 0.06; % vapor quality (mass fraction) 

%Main plant properties
MAIN_POWER = 1300; %MWe, power from main turbine of plant
THERMAL_POWER = 3500; %MWt, from steam generator
MDOTBASE = 1333; %kg/s, mass flow rate at main_power
MIN_LEVEL = 25; %percent, minimum turbine level as percent of main power
MIN_LOAD = (MIN_LEVEL/100)*MAIN_POWER; %MWe
p_topup = 1; %bar, low pressure makeup tank
h_topup = 206; %kJ/kg, low pressure makeup tank
sgh_output = 2770; %kJ/kg, outlet enthalpy at steam gen
sgh_input = sgh_output-(THERMAL_POWER*1000)/MDOTBASE; %kJ/kg, inlet enthalpy at steam gen
Number = 1; %number of cycles
POWER_REDUCTION = MAIN_POWER-MIN_LOAD;
MDOT_CHARGE = MDOTBASE-(MIN_LEVEL/100)*MDOTBASE; %kg/s, maximum charging mass flow rate to produce min_load at base case values

% %Sinusoidal price curve, amortization values, and other economic info
% life=40; %years, amortization period
% interest=0.07; %for amortization period
% period=6; %hours, price period
% peakAmplitude=25; %$/MWh
% avgElecPrice=34; %$/MWh
% start_cost = 11; %$/MW-start, from Neal
% var_OM = 8; %$/MWh, from Neal


%%runs steam_accumulator_separate code
%run either block 1 or block 2 (comment out the one you aren't using)

%%discharge cycle with net revenue and capital cost calculation
%block 1
discharge = 1; %discharge on
acc=steam_accumulator_separate_alt(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge);
acc.run_accumulator(POWER_SA, pow, discharge); %discharge
acc.charge(P0,X0,VTANK,MDOT_CHARGE,POWER_REDUCTION,POWER_SA,p_topup,h_topup,sgh_input,sgh_output);
%[netRevenue,CC,RC,RD,totalOM,totalCC]=acc.revenue(POWER_SA,ENERGY_SA,MAIN_POWER,MIN_LOAD,LTANK,life,interest,period,peakAmplitude,avgElecPrice,start_cost,var_OM);
   

%%store
%block 2
% discharge = 0; %discharge off (store)
% acc=steam_accumulator_separate_alt(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge);
% acc.run_accumulator(POWER_SA, pow, discharge); %store
    

%Insulation Thickness test:
% LossTotal = sum(acc.QLOSS); %kW lost over operation
% TimeOperated = acc.i*DT; % seconds
% LossPerHour = LossTotal*3600/TimeOperated; %kW/hour
% SystemPower =  POWER_SA * 1000 * 3600; %MW * kW/MW * sec/hour
% PercLoss = LossPerHour/SystemPower * 100; %Percent of Acc Power 
% disp(['Average % loss per hour is: ' num2str(PercLoss)])

LossTotal = sum(acc.QLOSS*(DT/3600)); %[kWh]
MaxEnergy = XSteam('hv_T',MaxTemp) * XSteam('rhoV_T',MaxTemp) * VTANK * 0.000277778; %[kWh]
PercLoss = LossTotal/MaxEnergy;
%disp(['Average % loss per hour is: ' num2str(PercLoss)])


%need to add more plots
acc.plots();
disp(acc.x(acc.N))
fprintf('Total run time = %.2f seconds.\n', toc(start_time));