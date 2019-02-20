T_HL = 275;
T_CL = 200;

SALT = 'Hitec';
%function [Th,Tc,Vh,Vc] = run_salt_thermo(T_HL,T_CL,SALT)
% This function out puts 4 arrays that depict the parameter value over
% time. The arrays are for hot and cold tanks temps and hot and cold tank
% volumes, respoectively. Inputs are only those that will be varied in the
% sensitivity analysis.
%==========================================================================================
% Inputs that are not changed in sensitivity analysis
%==========================================================================================
% Salt Tank Inputs  | Value             | Units       | Description/Notes
%-------------------|-------------------|-------------|-----------------------------
radius              = 6;                % m             both tanks have same dimensions
height              = 12;               % m             in Milestone 3 report a tank height of 12 was used
init_hot_vol_frac   = 0.75;             % -             total fluid volume equals the volume of one tank.  This is the fraction of the total fluid that's initially in the hot tank
P_in                = 50E6;             % J/s           rate at which THERMAL power is transferred in when charging
P_out               = 30E6;             % J/s           rate at which THERMAL power is transferred out when discharging
T_env               = 25;               % C             ambient temperature outside the tank

    C_p             = 1447;
    rho             = 1992;
    load HitecXLData
    T_data = HitecXLData(:,1);
    h_data = HitecXLData(:,7);



%Q_loss things

ksteel=41.8; %(W/mK)   
kinsul=0.079; %(W/mK)
Dtank = 2*radius; %m
Ttank = 0.2; %m
Tinsul = 0; %m
Tinf = 300; %(K)
D=Dtank+2*Tinsul; %total diameter (m)
h = 6.11; %30mph wind?
R1=log((Dtank/2)/((Dtank-2*Ttank)/2))/(2*pi()*ksteel); %conductive resistance of pipe (m*K/W)
R2=log((Dtank/2+Tinsul)/(Dtank/2))/(2*pi()*kinsul); %conductive resistance of insul (m*K/W)
Rconv=1/(h*pi()*D); %(m*K/W)
Rtotal=R1+R2+Rconv;



%-------------------------------------------------------------------------
% Initial Calculations
%-------------------------------------------------------------------------
volume              = pi*radius*radius*height;                % m3, tank volume
side_area           = 2*pi*radius*height + pi*radius^2;       % m2, area of tank sides + the top

%HXC is U, the total transfer coefficient
HXC = 1/(Rtotal*side_area);


%------------------------------------------------------------------------------------------
% Sets up events (discharge, nothing, charge)
%------------------------------------------------------------------------------------------
N_EVENTS                = 1;                %               number of events
events                  = [1]; % 0, 1 and 2 correspond to charge, store, and discharge, respectively
% can string together as many events 0, 1 and 2 as you want,
% but code will give garbage results if either tank becomes totally empty or if heat loss to environment is
% too high this can be fixed with some more coding but for now just keep an eye on the plots
endt_event              = [3600*100]; %...         % s             event 1 ends after 1800 s (30 minutes)
%     3600 ...         % s             event 2 ends at 3600 s (hence lasts 1800 s)
%     5400 ...         % s             sets the end time of each event
%     12000 ...
%     15000 ...
%     17600 ...
%     36000 ...
%     38400];
dt                      = 10;              % s             time step in seconds - I haven't tested convergence, could possibly jack this way up
tsteps_per_event(1)     = round(endt_event(1)/dt);       % number of time steps in the first event
endt_event(1)           = dt*tsteps_per_event(1);
nt                      = tsteps_per_event(1);
for i=2:N_EVENTS        % calculates total number of time steps for all events
    tsteps_per_event(i) = round((endt_event(i)-endt_event(i-1))/dt);% number of time steps in event i
    endt_event(i)       = endt_event(i-1) + dt*tsteps_per_event(i); % rounds entries here?
    nt                  = nt + tsteps_per_event(i);
end
t                       = (0:dt:(nt-1)*dt); % time array?

%==========================================================================================
% end of inputs
%==========================================================================================


%-------------------------------------------------------------------------
% Creates arrays of length nt for time, temperatures, volumes, and masses
%-------------------------------------------------------------------------
Th                      = [1,nt];
Tc                      = [1,nt];
Vh                      = [1,nt];
Vc                      = [1,nt];
Enthalpy_hot            = [1,nt];
Enthalpy_cold           = [1,nt];
mh                      = [1,nt];
mc                      = [1,nt];
Q_loss                  = [1,nt];
% mdotchg = zeros(nt,1)';
% mdotdis = zeros(nt,1)';


Th(1) = T_HL;                         % initial T in hot tank
Tc(1) = T_CL;                         % initial T in cold tank
Vh(1) = volume*init_hot_vol_frac;     % initial volume in hot tank
Vc(1) = volume*(1-init_hot_vol_frac); % initial volume in cold tank
Enthalpy_hot(1) = interp1(T_data,h_data,Th(1));
Enthalpy_cold(1) = interp1(T_data,h_data,Tc(1));
Enthalpy_hot_initial = Enthalpy_hot(1); 
Q_loss(1) = 0;
% mdotchg(1) = 0;
% mdotdis(1) = 0;
ix = 2;
% mt = [];

% calculated temperatures and volumes at the time steps for each case
% the formulas used for Th, Tc, Vh, and Vc can be found in the report -- derived by Schneider
for i=1:N_EVENTS
    if      events(i) == 0  % discharge
        for j=1:tsteps_per_event(i)
            Th(ix) = Th(ix-1)-dt*HXC*side_area*(Th(ix-1)-T_env)/(C_p*rho*Vh(ix-1));
            Tc(ix) = Tc(ix-1)-dt*HXC*side_area*(Tc(ix-1)-T_env)/(C_p*rho*Vc(ix-1));
            Vh(ix) = Vh(ix-1)-dt*P_out/(rho*C_p*(Th(ix-1)-T_CL));
            Vc(ix) = Vc(ix-1)+dt*P_out/(rho*C_p*(Th(ix-1)-T_CL));
            Enthalpy_hot(ix) = interp1(T_data,h_data,Th(ix));
            Enthalpy_cold(ix) = interp1(T_data,h_data,Tc(ix));

            % mt     = dt*(P_out/(rho*C_p*(Th(ix-1)-T_CL))); % Why is this the same as the second term in V eqn? Also this doesn't appear to be used later********
            ix     = ix+1;
        end % for int j
        ix = ix-2; % reset time step index after each event
    elseif   events(i) == 1  % just heat loss
        for j=1:tsteps_per_event(i)
            Th(ix) = Th(ix-1)-dt*HXC*side_area*(Th(ix-1)-T_env)/(C_p*rho*Vh(ix-1));
            Tc(ix) = Tc(ix-1)-dt*HXC*side_area*(Tc(ix-1)-T_env)/(C_p*rho*Vc(ix-1));
            Vh(ix) = Vh(ix-1);
            Vc(ix) = Vc(ix-1);
            Enthalpy_hot(ix) = interp1(T_data,h_data,Th(ix));
            Enthalpy_cold(ix) = interp1(T_data,h_data,Tc(ix));
            Tdelta=(Th(ix)+273)-Tinf; 
            Q_loss(ix) = Tdelta/Rtotal;
            ix     = ix+1;
            
        end %for int j
        ix = ix-2; % reset time step index after each event
    elseif  events(i) == 2  % charge
        for j=1:tsteps_per_event(i)
            Th(ix) = Th(ix-1) - dt*HXC*side_area*(Th(ix-1)-T_env)/(C_p*rho*Vh(ix-1)) + dt*P_in*(Th(ix-1)-T_HL)/(C_p*rho*Vh(ix-1)*(Tc(ix-1)-T_HL));
            Tc(ix) = Tc(ix-1) - dt*HXC*side_area*(Tc(ix-1)-T_env)/(C_p*rho*Vc(ix-1)) - dt*P_in*(Th(ix-1)-T_HL)/(C_p*rho*Vh(ix-1)*(Tc(ix-1)-T_HL));
            Vh(ix) = Vh(ix-1) - dt*P_in/(rho*C_p*(Tc(ix-1)-T_HL));
            Vc(ix) = Vc(ix-1) + dt*P_in/(rho*C_p*(Tc(ix-1)-T_HL));
            Enthalpy_hot(ix) = interp1(T_data,h_data,Th(ix));
            Enthalpy_cold(ix) = interp1(T_data,h_data,Tc(ix));
            ix     = ix+1;
        end %for int j
        ix = ix-2; % reset time step index after each event
    else
        break;
    end
end % for int i

MaxHotEnth = Enthalpy_hot(1);
MinHotEnth = Enthalpy_hot(ix+1);
TotalHeatLost = MaxHotEnth-MinHotEnth;
TimeOfProcess = nt*dt; %[s]
HoursOfProcess= TimeOfProcess/3600;
HeatLostperHour = TotalHeatLost/HoursOfProcess;
LossRate = HeatLostperHour/MaxHotEnth; % percent/hr
disp(LossRate)

%end
