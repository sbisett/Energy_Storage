%Calculates the maximum recoverable energy density of a steam accumulator
%Author:Seth Bisett

%Assume it is all steam at the presented maximum temperature, to find
%density. Next calculate the amount of pure steam in a 1 meter section of
%pipe (and the associated volume).Calculate the energy contained in that 
%steam mass, then apply the efficiency to find the maximum recoverable
%energy.

Tmax = 275;                                    % [C] Maximum Temperature of Accumulator
PipeOuterDiameter = 0.8128;                    % [m] 16 inch
PipeInnerDiameter = 0.8128-2*0.015875;         % [m] Outer Diameter - 2*Thickness
Efficiency = 0.33;                             % Reference plant efficiency

PipeInnerArea = PipeInnerDiameter^2*pi/4;      % [m^2]
PipeOuterArea = PipeOuterDiameter^2*pi/4;      % [m^2]
% VaporDensity = XSteam('rhoV_T',Tmax);          % [m^3/kg] From XSteam.m density of saturated vapor given temperature (Note:XSteam uses Celsius)
% PipeInnerVolume = PipeInnerArea*1;             % [m^3] Area times unit length
% PipeOuterVolume = PipeOuterArea*1;             % [m^3] Area times unit length
% VaporMass = VaporDensity*PipeInnerVolume;      % [kg]
% VaporEnthalpy = XSteam('hV_T',Tmax);           % [kJ\kg]
% VaporEnergy = VaporMass*VaporEnthalpy;         % [kJ Thermal]
% RecoverableEnergy_steam = VaporEnergy*Efficiency;    % [kJ Electric]
% RecoverablekWh_steam = 0.00027778*RecoverableEnergy_steam; % [kWh Electric]

LiquidDensity = XSteam('rhoL_T',Tmax);          % [m^3/kg] From XSteam.m density of saturated vapor given temperature (Note:XSteam uses Celsius)
PipeInnerVolume = PipeInnerArea*1;             % [m^3] Area times unit length
PipeOuterVolume = PipeOuterArea*1;             % [m^3] Area times unit length
LiquidMass = LiquidDensity*PipeInnerVolume;      % [kg]
LiquidEnthalpy = XSteam('hL_T',Tmax);           % [kJ\kg]
LiquidEnergy = LiquidMass*LiquidEnthalpy;         % [kJ Thermal]
RecoverableEnergy_steam = LiquidEnergy*Efficiency;    % [kJ Electric]
RecoverablekWh_steam = 0.00027778*RecoverableEnergy_steam; % [kWh Electric]



%disp(['The recoverable energy of the steam accumulator is ', num2str(RecoverablekWh_steam), ' [kWh/m] taking up a total volume of ', num2str(PipeOuterVolume), ' [m^3] including the pipe itself.'  ])


%This section for Salt 
Tmax = 275; 
Tmin = 200; 
%Get hitec data
load HitecXLData
T_data = HitecXLData(:,1);
h_data = HitecXLData(:,7);
rho_data = HitecXLData(:,3);
Enthalpy_high = interp1(T_data,h_data,Tmax); %[J/kg]
Enthalpy_low = interp1(T_data,h_data,Tmin);  %[J/kg]
Density = interp1(T_data,rho_data,Tmax); %[kg/m^3]
EnergyDensity_salt = (Enthalpy_high-Enthalpy_low)*Density/1000; % [kJ/m^3]
TankDiameter = 12; %[m]
TankHeight = 12; %[m]
TankVolume = TankHeight*pi*(TankDiameter^2)/4; %[m^3]
EnergyPerTank_salt = EnergyDensity_salt*TankVolume; %[kJ]
RecoverableEnergy_salt = EnergyPerTank_salt*Efficiency;
RecoverablekWh_salt =  0.00027778*RecoverableEnergy_salt; % [kWh Electric]



%disp(['The recoverable energy of the salt tank is ', num2str(RecoverablekWh_salt), ' [kWh] per tank taking up a total volume of ', num2str(TankVolume), ' [m^3].'  ])


%This section for 

ConcreteEnergyDensity_conc = 2.1e6; % [J/m^3]

%Cade is filling this out because the density per cube is variable