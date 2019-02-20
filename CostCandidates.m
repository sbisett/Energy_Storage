%This is the economic cost model for the project. It reads an input file with an .xls parser,and calculates scaled capital cost
%Everything else will ultimately be handled by the market model.
%Author: Seth Bisett, July 2018



%From EnergyDensityCalc.m
run('EnergyDensityCalc.m')
PipeEnergy = RecoverablekWh_steam*0.001; %[MWh]    %Calculated based on pipe length
TankEnergy = RecoverablekWh_salt*0.001; %[MWh]

%-------------------------------------------------------------------------%

%Setup to calculate scaling of cost data:
%TargetCost = CostScale(ScaleMethod,TargetSize,BaseSize,BaseCost,ScaleParameter)
%ScaleMethod:
    %1= Exponential , ScaleParameter = Scaling Exponent
    %2= Linear, ScaleParameter = 18
    %3= Integer, ScaleParameter = 1
    %4= Lin/Int Scaling, Scale Paramerer = # of tanks
    ExpScale = 1;
    LinScale = 2;
    IntScale = 3;
    LinIntScale = 4;
    
%automation loop    


%Inputs:
TechnologyUsed = 'salt'; %'steam' or 'salt' or 'concrete'
PowerBlockBase = 50; %[MW] Because the power block has a different reference size
InsulCost = 0.0002; %[MM $/m^3]
PipeLength = 1; %[m]
TankRadius = 6; %[m]
TankHeight = 2*TankRadius; %[m] Optimal Shape is h=2r for cylinder

CapitalCostReport = zeros(NumberOfCandidates,1);
for CandidateIteration = 1:NumberOfCandidates
    
    TargetPower = PowerArray(CandidateIteration);   %[MW]
    StorageHours = StorageTimeArray(CandidateIteration);
    TargetEnergy = TargetPower * StorageHours;  %[MWh] 
    InsulThickness = InsulationThicknessReport(CandidateIteration); %[m]
    
if strcmp(TechnologyUsed,'salt')
    CostData = xlsread('C:\Users\sethb\Desktop\Energy Storage\Model Data\Cost Data\SaltInputs.xlsx');
    %Tank number calculation
    NumTanks = ceil(TargetEnergy/TankEnergy);
    InsulArea = NumTanks*(2*pi*TankRadius*TankHeight+2*pi*TankRadius^2); %Insulation on all sides
   %Parameter----------Read/Scale----------------------------------------------------------------------Units------%
    BasePower =        CostData(1,1);                                                                % [MW]
    StorageTime =      CostData(2,1);                                                                % [Hours] 
    BaseEnergy =       BasePower * StorageTime;                                                      % [MWh]
    SaltCost =         CostScale(LinScale,TargetEnergy,BaseEnergy,CostData(4,1),1);                  % [MM $]
    TankCost =         CostScale(IntScale,TargetEnergy,BaseEnergy,CostData(5,1),1);                  % [MM $]
    InsulationCost =   InsulArea*InsulThickness*InsulCost;                                           % [MM $]
    FoundationCost =   CostScale(IntScale,TargetEnergy,BaseEnergy,CostData(7,1),1);                  % [MM $]
    HeatXCost =        CostScale(ExpScale,TargetPower,BasePower,CostData(8,1),CostData(7,2));        % [MM $]
    PumpCost =         CostScale(IntScale,TargetEnergy,BaseEnergy,CostData(10,1),CostData(10,2));    % [MM $]
    PowerBlockCost =   CostScale(ExpScale,TargetPower,PowerBlockBase,CostData(14,1),CostData(14,2)); % [MM $]
    PlantBalCost =     CostScale(ExpScale,TargetPower,BasePower,CostData(15,1),CostData(15,2));      % [MM $]
    GridConnCost =     CostScale(ExpScale,TargetPower,BasePower,CostData(16,1),CostData(16,2));      % [MM $]
    SysBalCost =       CostScale(ExpScale,TargetEnergy,BaseEnergy,CostData(11,1),CostData(11,2));    % [MM $]
    CapitalCosts =    [SaltCost,TankCost,FoundationCost,HeatXCost,PumpCost,PowerBlockCost,PlantBalCost,GridConnCost,SysBalCost,InsulationCost];
    
elseif strcmp(TechnologyUsed,'steam')
    CostData = xlsread('C:\Users\sethb\Desktop\Energy Storage\Model Data\Cost Data\SteamInputs.xlsx');
    %Pipe number calculation
    NumPipes = ceil( TargetEnergy/PipeEnergy);
    CR = NumPipes^(1/3);
    if CR == floor(CR)
       SizeUnits=CR^3;
    else
        CR = CR+1;
        SizeUnits = (CR)^3;
    end
    OuterBlocks = 6*CR^2;
    InsulArea = OuterBlocks * PipeLength;  %Should this be PipeLength or something else?
   %Parameter----------Read/Scale----------------------------------------------------------------------Units------%
    BasePower =        CostData(1,1);                                                                % [MW]
    StorageTime =      CostData(2,1);                                                                % [Hours] 
    BaseEnergy =       BasePower * StorageTime;                                                      % [MWh]
    PipeCost =         CostScale(IntScale,TargetEnergy,BaseEnergy,CostData(3,1),1);                  % [MM $]
    InsulationCost =   InsulArea*InsulThickness*InsulCost;                                           % [MM $] 
    FoundationCost =   CostScale(LinScale,TargetEnergy,BaseEnergy,CostData(5,1),1);                  % [MM $]
    SysBalCost =       CostScale(ExpScale,TargetEnergy,BaseEnergy,CostData(8,1),CostData(8,3));      % [MM $]
    PlantBalCost =     CostScale(ExpScale,TargetPower,BasePower,CostData(11,1),CostData(11,3));      % [MM $]
    PowerBlockCost =   CostScale(ExpScale,TargetPower,PowerBlockBase,CostData(10,1),CostData(10,3)); % [MM $]
    GridConnCost =     CostScale(ExpScale,TargetPower,BasePower,CostData(12,1),CostData(12,2));      % [MM $]  
    CapitalCosts =     [PipeCost,FoundationCost,SysBalCost,PlantBalCost,PowerBlockCost,GridConnCost,InsulationCost]; 
    
elseif strcmp(TechnologyUsed,'concrete')
    CostData = xlsread('C:\Users\sethb\Desktop\Energy Storage\Model Data\Cost Data\SteamInputs.xlsx'); %uses many of the same values from steam
    CostperVolume =    0.000285; %285 $/m^3
    CubeSideLength = 1;
    %Cube number calculation
    NumCubes = ceil( TargetEnergy/Cube0Energy);
    CR = NumCubes^(1/3);
    if CR == floor(CR)
       SizeUnits=CR^3;
    else
        CR = CR+1;
        SizeUnits = (CR)^3;
    end
    OuterBlocks = 6*CR^2;
    InsulArea = OuterBlocks * CubeSideLength;  
    
   %Parameter----------Read/Scale----------------------------------------------------------------------Units------%
    BasePower =        CostData(1,1);                                                                % [MW]
    StorageTime =      CostData(2,1);                                                                % [Hours] 
    BaseEnergy =       BasePower * StorageTime;                                                      % [MWh]
    PipeCost =         CostScale(IntScale,TargetEnergy,BaseEnergy,CostData(3,1),1);                  % [MM $]
    InsulationCost =   InsulArea*InsulThickness*InsulCost;                                           % [MM $] 
    FoundationCost =   CostScale(LinScale,TargetEnergy,BaseEnergy,CostData(5,1),1);                  % [MM $]
    SysBalCost =       CostScale(ExpScale,TargetEnergy,BaseEnergy,CostData(8,1),CostData(8,3));      % [MM $]
    PlantBalCost =     CostScale(ExpScale,TargetPower,BasePower,CostData(11,1),CostData(11,3));      % [MM $]
    PowerBlockCost =   CostScale(ExpScale,TargetPower,PowerBlockBase,CostData(10,1),CostData(10,3)); % [MM $]
    GridConnCost =     CostScale(ExpScale,TargetPower,BasePower,CostData(12,1),CostData(12,2));      % [MM $]  
    ConcreteCost =     Units*VolumeperUnit*CostperVolume;                                            % [MM $]
    CapitalCosts =     [PipeCost,FoundationCost,SysBalCost,PlantBalCost,PowerBlockCost,GridConnCost,InsulationCost,ConcreteCost]; 
    
else 
    error('Technology not correctly specified!')
end

TotalCapitalCost = sum(CapitalCosts);
TotalCost_CC_OM = TotalCapitalCost*(4/3);
disp(['The total capital cost of this ' TechnologyUsed ' system is: ' num2str(TotalCost_CC_OM) ' million dollars.'])
CapitalCostReport(CandidateIteration) = 1000000*TotalCost_CC_OM;

end

save ('CapitalCostReport_salt','CapitalCostReport')