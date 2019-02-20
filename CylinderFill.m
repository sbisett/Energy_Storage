%Tries to calculate how to best fit pipes in a cylinder configuration
%Author: Seth Bisett

%Take energy density from EnergyDensityCalc.m and come up with how many
%pipes for a 100MWh accumulator.

RecoverablekWh = 4.3887; %[kWh/m]
Target = 100000; %[kWh]

PipeOuterDiameter = 0.8128;                    % [m] 16 inch
PipeOuterArea = PipeOuterDiameter^2*pi/4;      % [m^2]
PipeOuterVolume = PipeOuterArea*1;             % [m^3] Area times unit length

Length = Target/RecoverablekWh;

%Best surface area to volume for a pipe is h = 2r
%Rough assumption, pipes fit perfectly area is PipeOuterArea with no
%extraneous
TotalVolume = Length*PipeOuterArea;
%Volume of a cylinder = h*pi*r^2, h*pi*h^2/4 for 'ideal' cylinder
hIdeal = (TotalVolume*4/pi)^(1/3);
rIdeal = (TotalVolume/(hIdeal*pi))^0.5;

disp(['Our ',num2str(Target/1000) ' [MWh] accumulator will take up *at least* ', num2str(TotalVolume) ' [m^3] and have length and radius of ', num2str([hIdeal,rIdeal]),' [m] under impossible packing conditions.'])
