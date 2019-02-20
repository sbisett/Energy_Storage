clear variables
close all
%run this script first
start_time = tic;



%Properties for runs:


%SA properties

POWER_SA = 1400; %MW, maximum electric power provided by steam accumulator
T_STORE = 24*3600; %s, storage time
LossRateTarget = 0.1; %[%]
InsulationThickness = 0; %[m]


%accumulator initial thermo properties
P0 = 70; % bar, initial pressure
X0 = 0.06; % vapor quality (mass fraction) 


%set stuff up ahead of time:
%SA properties
T_MAX = 3; %hr, time at maximum SA power
D_RAMP_RATE = 1.67; %percent/min, discharge ramp rate (SA turbine, % of max power/min)
POWER_SA_INITIAL = 0; %MW, cold start
T_RAMP = 60*((POWER_SA-POWER_SA_INITIAL)/((D_RAMP_RATE/100)*POWER_SA)); %s, time to ramp up to max power
T_END = T_RAMP+T_MAX*3600; %s, discharge time including ramp up and time at max power
ENERGY_SA = 0.5*T_RAMP*(POWER_SA-POWER_SA_INITIAL)/3600+POWER_SA*T_MAX; %MWh, electric energy provided by SA in one discharge cycle
DT = 30;  % s, time step
RTANK = 0.4064; % m (16 inches)
MaxTemp = 285; %Used in loss rate/ hr calculation

t_length = floor(T_RAMP/DT); %length of ramping power vector
pow = zeros(t_length,1);
for t = 1:t_length
    pow(t) = t*DT*((D_RAMP_RATE/100)*POWER_SA)/60; %MW, power as SA turbine is ramping up
end

%Automation of runs setup:

Candidates = xlsread('C:\Users\sethb\Desktop\Energy Storage\Candidates\Heat storage candidates matrix STEAM LowLoss');
LossRateTargetArray = Candidates(:,1);ans
PowerArray = Candidates(:,2);
StorageTimeArray = Candidates(:,3);
NumberOfCandidates = length(PowerArray);
CapitalCostOutputs = zeros(NumberOfCandidates,1);
InsulationThicknessReport = zeros(NumberOfCandidates,1);
RealLossRateReport = zeros(NumberOfCandidates,1);
%Automation loop
for CandidateIteration = 1:NumberOfCandidates
    
    %set variables
    T_STORE = StorageTimeArray(CandidateIteration)*3600; %s, storage time
    if T_STORE == 4*3600
        LTANK = 450000;
    elseif T_STORE == 8*3600
        LTANK = 880000;
    elseif T_STORE == 12*3600
        LTANK = 1380000;
    elseif T_STORE == 16*3600
        LTANK = 1900000;
    elseif T_STORE == 20*3600
        LTANK = 2500000;
    elseif T_STORE == 24*3600
        LTANK = 3150000;
    end
    POWER_SA = PowerArray(CandidateIteration); %MW, maximum electric power provided by steam accumulator
    LossRateTarget = LossRateTargetArray(CandidateIteration); %[%]
    InsulationThickness = 0; %[m]
    VTANK = LTANK.*pi.*RTANK^2.;  % m3
    
    %run simulation
    discharge = 0; %discharge off (store)
    acc=steam_accumulator_separate(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge,InsulationThickness);
    acc.run_accumulator(POWER_SA, pow, discharge,InsulationThickness); %store
    LossTotal = sum(acc.QLOSS*(DT/3600)); %[kWh]
    MaxEnergy = XSteam('hv_T',MaxTemp) * XSteam('rhoV_T',MaxTemp) * VTANK * 0.000277778; %[kWh]
    PercLossperHour = LossTotal/MaxEnergy/StorageTimeArray(CandidateIteration)*100;
    disp(['Candidate: ', num2str(CandidateIteration)])
    disp(['Time ',num2str(StorageTimeArray(CandidateIteration)) , ' Length ',num2str(LTANK) ])
    disp(['Insulation Thickness:',num2str(InsulationThickness)])
    disp(['%loss/hr ',num2str(PercLossperHour)])
    while (PercLossperHour - LossRateTarget > 0)
        InsulationThickness = InsulationThickness + 0.5;
        %Insulation Thickness test:
        LossTotal = sum(acc.QLOSS*(DT/3600)); %[kWh]
        MaxEnergy = XSteam('hv_T',MaxTemp) * XSteam('rhoV_T',MaxTemp) * VTANK * 0.000277778; %[kWh]
        PercLossperHour = LossTotal/MaxEnergy/StorageTimeArray(CandidateIteration)*100;
        acc=steam_accumulator_separate(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge,InsulationThickness);
        acc.run_accumulator(POWER_SA, pow, discharge,InsulationThickness); %store
        disp(['Candidate: ', num2str(CandidateIteration)])
        disp(['Time ',num2str(StorageTimeArray(CandidateIteration)) , ' Length ',num2str(LTANK) ])
        disp(['Insulation Thickness: ',num2str(InsulationThickness)])
        disp(['%loss/hr: ',num2str(PercLossperHour)])
    end
    InsulationThicknessReport(CandidateIteration) = InsulationThickness;
    RealLossRateReport(CandidateIteration) = PercLossperHour;
end

save('SteamCandidateRunData_partial')
fprintf('Total run time = %.2f seconds.\n', toc(start_time));