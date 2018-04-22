classdef steam_accumulator_separate < handle
    properties
        p = []; %bar
        m1 = []; %kg, liquid
        m2 = []; %kg, vapor
        Vol1 = []; %m^3
        Vol2 = []; %m^3
        v1 = []; %m^3/kg
        v2= []; %m^3/kg
        Vol_total; %m^3
        rho1 = []; %kg/m^3
        rho2 = []; %kg/m^3
        mix_rho = [];
        h1 = []; %kJ/kg
        h2 = []; %kJ/kg
        h2in = []; %kJ/kg
        h2out = []; %kJ/kg
        x = []; %quality, (mass fracton of vapor)
        t1 = []; %C
        t2 = []; %C
        m1b = []; %kg/s
        m2b = []; %kg/s
        mh1b = []; %kJ/s
        mh2b = [];
        mpt1 = [];   % kg/s
        mpt2 = [];   % kg/s; mpt2=-mpt1
        m_total; %kg
        m_original; %kg
        QLOSS = []; %kW
        qloss1 = []; %kW
        qloss2 = []; %kW
        VAP_IN; %kg/s
        VAP_OUT; %kg/s
        LIQ_IN; %kg/s
        LIQ_OUT; %kg/s
        PVAP_IN; %bar
        PLIQ_IN; %bar
        HVAP_IN; %kJ/kg
        HLIQ_IN; %kJ/kg
        N; %array length, used for initial setup
        time = []; %time array
        i; % counter
        time_step; %s
        discharge_time; %s
        store_time; %s
        ramp_time; %s
        length; %m, pipe length
        radius; %m, pipe radius
        Dh; %kJ
        Dp; %bar
        dv1dh;
        dv2dh;
        dv1dp;
        dv2dp;
        WTURB;
        r;
        me;
        mc;
        hL_sat; %kJ/kg
        hV_sat; %kJ/kg
        Area_int;
        Qint; %kW, interphase heat transfer
        charge_time; %s
        efficiency_r;
        thermo_eff; %efficiency of ideal rankine cycle connected to SA turbine
        P_disch_final; %bar, pressure at end of discharge
        X_disch_final; %quality at end of discharge
        POWER_LOSS_AVERAGE_DISCHARGE; %kW, due to heat loss
        POWER_LOSS_AVERAGE_STORE; %kW, due to heat loss
        PERCENT_LOSS_STORE; %percent, loss due to heat loss during storage
        ENERGY_LOSS_DISCHARGE; %kJ
        ENERGY_LOSS_STORE; %kJ
        PERCENT_LOSS_DISCHARGE; %percent, loss due to heat loss during discharge.
        HH_initial; %kJ, initial enthalpy of fully charged steam accumulator
        charge_ramp; %percent/s, charging ramp rate
    end
    properties(Constant)
       PASCALS_PER_BAR = 1.e5;
       EPSILON = .01;
       TAU = 85; %s
       deltaH=-.01; %kJ/kcg, small change in enthalpy used to evaluate numerical derivative
       deltaP=-.01; %bar, small change in pressure used to evaluate numerical derivative
       Coeff_int= 910; %Interphase heat transfer coefficient[W/m^2K]
    end
    methods
        function object = steam_accumulator_separate(T_DISCHARGE,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge)
            object.i=1;
            object.time_step=DT;
            object.discharge_time=T_DISCHARGE;
            object.ramp_time=T_RAMP;
            object.store_time=T_STORE;
            object.length = LTANK;
            object.radius = RTANK;
            object.setup_arrays(discharge);
            object.initial_conditions(X0, P0, VTANK);
        end
        function run_accumulator(object, POWER_SA, pow, discharge)
            object.VAP_IN = 0; %during discharge 
            object.LIQ_IN = 0;
            object.LIQ_OUT = 0;
            object.PVAP_IN = 70.; %bar (not utilized during discharge)
            object.PLIQ_IN = 70.; %bar (not utilized during discharge)
            object.HVAP_IN = XSteam('hV_p',object.PVAP_IN); %kj/kg
            object.HLIQ_IN = XSteam('hL_p',object.PLIQ_IN); %kj/kg
            
            if discharge == 1
            [object.VAP_OUT(object.i), object.thermo_eff(object.i)]=rankine(object, pow(object.i), object.p(object.i));
            else
                object.VAP_OUT(object.i)=0; %sets mass flow rate to zero during storage
            end
            for s=1:(object.N-1)
                object.m1b(object.i)=object.LIQ_IN-object.LIQ_OUT; %kg/s
                object.m2b(object.i)=object.VAP_IN-object.VAP_OUT(object.i); %kg/s
                object.mh1b(object.i)=object.LIQ_IN*object.HLIQ_IN-object.LIQ_OUT*object.h1(object.i); %kW
                object.mh2b(object.i)=object.VAP_IN*object.HVAP_IN-object.VAP_OUT(object.i)*object.h2(object.i); %kW

                object.hL_sat = XSteam('hL_p',object.p(object.i)); %saturated liquid enthalpy (kJ/kg)
                object.hV_sat = XSteam('hV_p',object.p(object.i)); %saturated vapor enthalpy (kJ/kg)

                object.r = object.hV_sat - object.hL_sat;     % latent heat of vaporization
                if(object.h1(object.i) > object.hL_sat)
                    object.mc = 0.;
                    object.me(object.i) = (object.rho1(object.i)*object.Vol1(object.i)*(object.h1(object.i)-object.hL_sat))/(object.TAU*object.r);
                else
                    object.me(object.i) = 0.;
                    object.mc = object.rho1(object.i)*object.Vol1(object.i)*(object.hL_sat-object.h1(object.i))/(object.TAU*object.r);
                end

                object.Area_int(object.i+1) = 2*object.radius*cos(asin(((object.Vol1(object.i)/object.Vol_total(object.i))-object.radius)/object.radius))*object.length; %Minimum Interface Area between the two phases[m^2/m^3]
                object.Qint(object.i+1) =((object.Coeff_int*object.Area_int(object.i+1))*(object.t2(object.i)-object.t1(object.i))*object.Vol1(object.i))/1000; %kW, Rate of interphase heat transfer from vapor to liquid [kW]
                %object.Qint(object.i+1) = 0; %kW
                
                object.mpt1(object.i)=object.mc-object.me(object.i); %kg/s, liquid mass change due to evaporation and condensation
                object.mpt2(object.i)=object.me(object.i)-object.mc; %kg/s, vapor mass change due to evaporation and condensation
                object.m1(object.i+1)=object.m1(object.i)+(object.m1b(object.i)+object.mpt1(object.i))*object.time_step; %liquid mass balance (kg)
                object.m2(object.i+1)=object.m2(object.i)+(object.m2b(object.i)+object.mpt2(object.i))*object.time_step; %vapor mass balance (kg)

                object.dv1dh(object.i)=(XSteam('v_ph',object.p(object.i),object.h1(object.i)+object.Dh)-XSteam('v_ph',object.p(object.i),object.h1(object.i)))/object.Dh; % (m3/kJ), change in v1 per change in h at constant p
                object.dv2dh(object.i)=(XSteam('v_ph',object.p(object.i),object.h2(object.i)+object.Dh)-XSteam('v_ph',object.p(object.i),object.h2(object.i)))/object.Dh; % (m3/kJ)
                object.dv1dp(object.i)=(XSteam('v_ph',object.p(object.i)+object.Dp,object.h1(object.i))-XSteam('v_ph',object.p(object.i),object.h1(object.i)))/(object.Dp*object.PASCALS_PER_BAR); %(m5/(N*kg)) change in v1 per change in p at constant h
                object.dv2dp(object.i)=(XSteam('v_ph',object.p(object.i)+object.Dp,object.h2(object.i))-XSteam('v_ph',object.p(object.i),object.h2(object.i)))/(object.Dp*object.PASCALS_PER_BAR);

                term1 = (object.h1(object.i)*object.dv1dh(object.i)-object.v1(object.i))*(object.m1(object.i+1)-object.m1(object.i))/object.time_step; %m3/s
                term2 = (object.h2(object.i)*object.dv2dh(object.i)-object.v2(object.i))*(object.m2(object.i+1)-object.m2(object.i))/object.time_step; %m3/s
                term3 = object.dv1dh(object.i)*(object.mh1b(object.i)+object.mpt1(object.i)*object.hV_sat-object.qloss1(object.i)+object.Qint(object.i)); %m3/s
                term4 = object.dv2dh(object.i)*(object.mh2b(object.i)+object.mpt2(object.i)*object.hV_sat-object.qloss2(object.i)-object.Qint(object.i)); %m3/s
                term5 = object.m1(object.i)*(object.dv1dp(object.i)+object.v1(object.i)*object.dv1dh(object.i)*.001); %m5/N
                term6 = object.m2(object.i)*(object.dv2dp(object.i)+object.v2(object.i)*object.dv2dh(object.i)*.001); %m5/N

                %% final values at end of disch are here
                object.p(object.i+1) = ((object.time_step*((term1+term2-term3-term4)/(term5+term6)))/object.PASCALS_PER_BAR)+object.p(object.i); %bar
                object.h1(object.i+1) = (object.mh1b(object.i)+object.mpt1(object.i)*object.hV_sat+(object.Vol1(object.i)/object.time_step)*object.PASCALS_PER_BAR*(object.p(object.i+1)-object.p(object.i))*.001-(object.h1(object.i)/object.time_step)*(object.m1(object.i+1)-object.m1(object.i)))*(object.time_step/object.m1(object.i))+object.h1(object.i)+object.Qint(object.i)*(object.time_step/object.m1(object.i)); %kJ/kg
                object.h2(object.i+1) = (object.mh2b(object.i)+object.mpt2(object.i)*object.hV_sat+(object.Vol2(object.i)/object.time_step)*object.PASCALS_PER_BAR*(object.p(object.i+1)-object.p(object.i))*.001-(object.h2(object.i)/object.time_step)*(object.m2(object.i+1)-object.m2(object.i)))*(object.time_step/object.m2(object.i))+object.h2(object.i)-object.Qint(object.i)*(object.time_step/object.m1(object.i)); %kJ/kg
                object.t1(object.i+1) = XSteam('T_ph',object.p(object.i+1),object.h1(object.i+1));
                object.t2(object.i+1) = XSteam('T_ph',object.p(object.i+1),object.h2(object.i+1));
                object.rho1(object.i+1) = XSteam('rho_ph',object.p(object.i+1),object.h1(object.i+1));
                object.rho2(object.i+1) = XSteam('rho_ph',object.p(object.i+1),object.h2(object.i+1));
                object.Vol1(object.i+1) = object.m1(object.i+1)/object.rho1(object.i+1);
                object.Vol2(object.i+1) = object.m2(object.i+1)/object.rho2(object.i+1); 
                object.Vol_total(object.i+1)=object.Vol1(object.i+1)+object.Vol2(object.i+1); 
                object.v1(object.i+1) = 1/object.rho1(object.i+1);
                object.v2(object.i+1) = 1/object.rho2(object.i+1);
                object.x(object.i+1) = object.m2(object.i+1)/(object.m1(object.i+1)+object.m2(object.i+1));
                object.QLOSS(object.i+1) = object.length*Q_loss(object,(object.t1(object.i+1)+object.t2(object.i+1))/2); %kW
                object.qloss1(object.i+1)=(object.Vol1(object.i+1)/(object.Vol1(object.i+1)+object.Vol2(object.i+1)))*object.QLOSS(object.i+1); %(kW)
                object.qloss2(object.i+1)=(object.Vol2(object.i+1)/(object.Vol1(object.i+1)+object.Vol2(object.i+1)))*object.QLOSS(object.i+1); %(kW)
                
                if isnan(object.p(object.i+1))
                    disp('Error. Increase pipe length');
                    break
                end
                
                if discharge == 1 %if discharging
                    if object.i<=floor(object.ramp_time/object.time_step) %if turbine is ramping
                    [object.VAP_OUT(object.i+1), object.thermo_eff(object.i+1)]=rankine(object, pow(object.i), object.p(object.i+1));          
                    else %else, turbine is at max power
                    [object.VAP_OUT(object.i+1), object.thermo_eff(object.i+1)]=rankine(object, POWER_SA, object.p(object.i+1));     
                    end
                else %else storing
                    object.VAP_OUT(object.i+1)=0;
                end
                object.i = object.i+1;
            end
            sum=0;
            for count=1:length(object.QLOSS)
            sum = sum + object.QLOSS(count)*object.time_step;
            end
            if discharge == 1
                object.POWER_LOSS_AVERAGE_DISCHARGE = sum*(1/(object.discharge_time)); %kW
                object.ENERGY_LOSS_DISCHARGE = sum; %kJ
                object.PERCENT_LOSS_DISCHARGE = (object.ENERGY_LOSS_DISCHARGE/object.HH_initial)*100; %percent loss during discharge due to heat loss
            else
                object.POWER_LOSS_AVERAGE_STORE = sum*(1/object.store_time); %kW
                object.ENERGY_LOSS_STORE = sum; %kJ
                object.PERCENT_LOSS_STORE = (object.ENERGY_LOSS_STORE/object.HH_initial)*100 %percent loss during storage due to heat loss
            end
            object.P_disch_final = object.p(object.i-1); %bar, final value after discharge
            object.X_disch_final = object.x(object.i-1); %final vapor quality
        
        end
        %%%%%charging block
        function object = charge(object,P0,X0,VTANK,MDOT_CHARGE,POWER_REDUCTION,POWER_SA,p_topup,h_topup,sgh_input,sgh_output)
           
            P_INITIAL = object.P_disch_final; % bar, final value after discharge
            P_END = P0; %bar
            m_target = object.m_total;  % mass of full accumulator

            X_INITIAL = object.X_disch_final; % initial vapor quality (mass fraction) 
            X_END = X0; % final target vapor quality
            rho1_INITIAL=XSteam('rhoL_p',P_INITIAL); %kg/m3
            rho2_INITIAL=XSteam('rhoV_p',P_INITIAL); %kg/m3
            t1_INITIAL=XSteam('Tsat_p',P_INITIAL);
            t2_INITIAL=t1_INITIAL;

            mixture_rho=1./(X_INITIAL/rho2_INITIAL+(1.-X_INITIAL)/rho1_INITIAL);
            mass_total=mixture_rho*VTANK; %total mass in tank

            m1_INITIAL = mass_total*(1-X_INITIAL); %total liquid mass in tank (kg)
            m2_INITIAL = mass_total*X_INITIAL; %total vapor mass in tank (kg)

            v1_INITIAL=1/rho1_INITIAL; %liquid specific volume in tank (m3/kg)
            v2_INITIAL=1/rho2_INITIAL; %vapor specific volume in tank (m3/kg)
            h1_INITIAL = XSteam('h_pT',P_INITIAL,t1_INITIAL-object.EPSILON); %(kj/kg)
            h2_INITIAL = XSteam('h_pT',P_INITIAL,t2_INITIAL+object.EPSILON);


            %% mix masses to find accum starting point, taking makeup mass from LP turbine discharge

            HH_accum_initial = m1_INITIAL*h1_INITIAL+m2_INITIAL*h2_INITIAL; % kJ, that's just for the residue in the accumulator, next add top-up
            m_topup = m_target - (m1_INITIAL+m2_INITIAL); % kg, amount of water needed to return to target mass
            v_topup = XSteam('v_ph',p_topup,h_topup); %m^3/kg, specific volume of top-up water

            HH_topup_initial = h_topup*m_topup; % kJ, extensive enthalpy of top up water
            HH_inital = HH_topup_initial+HH_accum_initial; % kJ, total extensive enthalpy of accum + top up water at start of recharge
            t_final = XSteam('Tsat_p',P_END); %C, final temp after charging
            h1_final = XSteam('hL_p',P_END); %kj/kg
            h2_final = XSteam('hV_p',P_END); %kJ/kg
            v1_final = XSteam('vL_p',P_END); %m^3/kg
            v2_final = XSteam('vV_p',P_END); %m^3/kg

            HH1_final = h1_final*m_target*(1-X_END); %kJ
            HH2_final = h2_final*m_target*X_END; %kJ
            HH_final = HH1_final+HH2_final; % kJ, total enthalpy of fully recharged accumulator(our target)
            HH_needed = HH_final - HH_inital; % kJ, this is how much enthalpy neeeds to be provided.  It will lead us to the charging time

            % other component: work done to pressurize accum to P0 

            pV_initial = (m1_INITIAL+m2_INITIAL)*P_INITIAL*object.PASCALS_PER_BAR*(v1_INITIAL*(1-X_INITIAL)+v2_INITIAL*X_INITIAL)/1000+p_topup*v_topup*m_topup*object.PASCALS_PER_BAR/1000; %kJ
            pV_final = P_END*object.PASCALS_PER_BAR*m_target*(v1_final*(1-X_END)+v2_final*X_END)/1000; %kJ

            pV_delta = pV_final - pV_initial; %kJ
            Total_input = HH_needed + pV_delta; %kJ

            object.charge_time = Total_input/MDOT_CHARGE/(sgh_output-sgh_input); %s,  total input [kJ] required from main steam, mDOTCHARGE a fixed value
            %sgh_output and sgh_input are inlet and outlet conditions of steam
            %generator
            object.charge_ramp=((sgh_output*MDOT_CHARGE)/object.HH_initial)*100; %percent/s, charge ramp rate: enthalpy added as percent of fully charged enthalpy
            %efficiency ratio calc. Check these. Charge time should be >
            %discharge time
            charge_energy_loss = object.charge_time*POWER_REDUCTION;  % dynamically grab power reduction during chg.
            discharge_energy_boost = object.discharge_time*POWER_SA;  %MJ, disch. power set by user
            object.efficiency_r = discharge_energy_boost/charge_energy_loss;
            
        end
        
        %Kayla's section 
        function [netRevenue,CC,RC,RD,totalOM,totalCC]=revenue(object,POWER_SA,ENERGY_SA,MAIN_POWER,MIN_LOAD,LTANK,life,interest,period,peakAmplitude,avgElecPrice,caseNumber, hotCyclesPerYear, warmCyclesPerYear, coldCyclesPerYear, var_om)   
           
            %this section mainly determines the charge and discharge time
            %discharge time is calculated a value input in RUN_SA_SEPARATE
            %charge time for the SA is calculated with line 276, and charge time for the molten salt is calculated with line 277
            %comment/uncomment one of the charge time lines depending on which model is being used
            %length must be connected to capital cost through pipe and insulation costs
            Y=(1/period)*24*365; %storage cycles per year
            d_t = object.discharge_time/3600; %hours, discharge time
            %c_t = object.charge_time/3600 %hours, charge time for SA
            c_t = (POWER_SA + 7.61232)/(POWER_SA-7.61232)*d_t %hours, charge time for salt
             
            %case number input in RUN_SA_SEPARATE assinged to local variable. Text file name also assigned
            %SA text file: SA_cost_input.txt
            %salt text file: salt_cost_input.txt
            %change fileName as indicated below
            caseName = 'CASE' + string(caseNumber);
            fileName = 'salt_cost_input.txt'; %change file name in otder to test different ones 
            
            %this section searches for the case name in each line of the text file
            %when the case name is found in a line, it assignes all the following text file lines in that case to a MATLAB variable 
            fileID = fopen(fileName);
                while ~feof(fileID)
                    tline = fgetl(fileID);
                    if contains(tline, caseName) 
                        caseLine = tline;
                        powerEnergyLine = fgetl(fileID);
                        costLine = fgetl(fileID);
                        scaleFactorLine = fgetl(fileID);
                        startCostLine = fgetl(fileID);
                    end
                end
                fclose(fileID);
            
           
            %in the following sections (lines 307-328), the lines assigned above are parsed to find numerical values of interest (power/energy costs, scale factors, etc)
            %numerical values are put into an array 
            %NOTE: the spaces before the units (like in line 307) are exact (based on the spaces in the text file -- don't change)
                
            %dive power energy line into array of doubles
            powerEnergyLine2 = extractBefore(powerEnergyLine, '     % MW MWh');
            powerEnergySpaces = regexp(powerEnergyLine2, ' ');
            powerEnergyStringArray = [extractBefore(powerEnergyLine2, powerEnergySpaces(1)) extractAfter(powerEnergyLine2, powerEnergySpaces(1))];
            powerEnergyDoubleArray = str2double(powerEnergyStringArray);   
            
            %divide cost line into array of doubles
            costLine2 = extractBefore(costLine,'   MM$ MM$ MM$/yr MM$/yr');
            costSpaces = regexp(costLine2, ' ');
            costStringArray = [extractBefore(costLine2, costSpaces(1)) extractBetween(costLine2, costSpaces(1)+1, costSpaces(2)-1) extractBetween(costLine2, costSpaces(2)+1, costSpaces(3)-1) extractAfter(costLine2, costSpaces(3))];
            costDoubleArray = str2double(costStringArray);    
            
            %divide scale factor line into array of doubles
            scaleFactorLine2 = extractBefore(scaleFactorLine,'		    % % % %');   
            scaleFactorSpaces = regexp(scaleFactorLine2, ' ');
            scaleFactorStringArray = [extractBefore(scaleFactorLine2, scaleFactorSpaces(1)) extractBetween(scaleFactorLine2, scaleFactorSpaces(1)+1, scaleFactorSpaces(2)-1) extractBetween(scaleFactorLine2, scaleFactorSpaces(2)+1, scaleFactorSpaces(3)-1) extractAfter(scaleFactorLine2, scaleFactorSpaces(3))];
            scaleFactorDoubleArray = str2double(scaleFactorStringArray);
            
            %divide start cost line into array of doubles
            startCostLine2 = extractBefore(startCostLine, '      $/MW-Cycle');
            startCostSpaces = regexp(startCostLine2, ' ');
            startCostStringArray = [extractBefore(startCostLine2, startCostSpaces(1)) extractBetween(startCostLine2, startCostSpaces(1)+1, startCostSpaces(2)-1) extractAfter(startCostLine2, startCostSpaces(2))];
            startCostDoubleArray = str2double(startCostStringArray);
            
            
            %in these sections (lines 334-358), values from the above arrays are assigned to individual variables
       
            %reference power and energy assigned
            Pref = powerEnergyDoubleArray(1); %reference power
            Eref = powerEnergyDoubleArray(2); %reference energy
            
            %power/energy cost and OM assigned 
            Cp = costDoubleArray(1); %cost (power)
            Ce = costDoubleArray(2); %cost (energy)
            Op = costDoubleArray(3); %OM (power)
            Oe = costDoubleArray(4); %OM (energy)
            
            %scale factors assigned for power/energy cost and OM 
            n_Cp = scaleFactorDoubleArray(1); %scaleFactor cost (power)
            n_Ce = scaleFactorDoubleArray(2); %scaleFactor cost (energy)
            n_Op = scaleFactorDoubleArray(3); %scaleFactor OM (power)
            n_Oe = scaleFactorDoubleArray(4); %scaleFactor cost (energy)
            
            %Cp Ce Op Oe scaled and assinged 
            Cp_scaled = Cp * (POWER_SA/Pref)^n_Cp; %cost scaled (power)
            Ce_scaled = Ce * (ENERGY_SA/Eref)^n_Ce; %cost scaled (energy)
            Op_scaled = Op * (POWER_SA/Pref)^n_Op; %OM scaled (power)
            Oe_scaled = Oe * (ENERGY_SA/Eref)^n_Oe; %OM sclaed (energy)
         
            %cycles OM calcaulated 
            coldStart = startCostDoubleArray(1); %$/MW-Cycle
            warmStart = startCostDoubleArray(2); %$/MW-Cycle
            hotStart = startCostDoubleArray(3); %$/MW-Cycle
            cyclingCost = ((coldStart * coldCyclesPerYear + warmStart * warmCyclesPerYear + hotStart * hotCyclesPerYear) * POWER_SA)/1000000; %MM$/year
            
            %calculations to determine ADC, ACP, DP
            c1=(3/4)*period-c_t/2; %hr, charge time integral lower bound
            c2=(3/4)*period+c_t/2; %hr, charge time integral upper bound
            d1=(period/4)-d_t/2; %hr, discharge time integral lower bound
            d2=(period/4)+d_t/2; %hr, discharge time integral upper bound
            y=@(t)peakAmplitude*sin((2*pi()*t)/period)+avgElecPrice; 
            intC=integral(y,c1,c2); %$/MW
            intD=integral(y,d1,d2); %$/MW
            ADP= intD/(d2-d1); %$/MWh, Average discharge price
            ACP= intC/(c2-c1); %$/MWh, Average charge price
            DP=ADP-ACP ; %$/MWh, delta price
  
            %final costs/revenues caluclated and displayed 
            %These need to be updated to account for ramping.
            RC=ACP*c_t*Y*(MAIN_POWER-MIN_LOAD)/10^6; %MM$/year, forgone revenue from charging
            RD=ADP*d_t*Y*POWER_SA/10^6; %MM$/year, revenue from discharging
            totalCC=Cp_scaled+Ce_scaled %Million $, total overnight capital cost
                disp('   MM$')
            CC=totalCC*(interest+(interest/((1+interest)^life-1))); %MM$/year, amortized capital cost
            fixed_om = (Oe_scaled+Op_scaled)*(10^6)*(1/POWER_SA)*(1/1000) %$/kw-year
                disp('   $/kw-year')
            totalOM = Op_scaled + Oe_scaled + cyclingCost + var_om *ENERGY_SA*Y/1000000 %MM$/year
                disp('   MM$/year')
            startCost = cyclingCost*1000000/(POWER_SA*(hotCyclesPerYear+warmCyclesPerYear+coldCyclesPerYear)) %$/MW-start
                disp('   $/MW-start')
            netRevenue=RD-RC-CC-totalOM %MM$/year, revenue provided by the addition of the accumulator
                disp('   MM$/year')
                disp(' ')
        end
        
        function object = plots(object)
            figure (1)
            plot(object.time,object.QLOSS)
            xlabel('time [s]')
            ylabel('Heat loss [kW]')
            
            figure (2)
            plot(object.time,object.p)
            xlabel('time [s]')
            ylabel('Pressure [bar]')
            
            figure (3)
            plot(object.time,object.x)
            xlabel('time [s]')
            ylabel('quality')
        end
      function Q_loss = Q_loss(~,internal_temp)
            %this function takes in the temperature inside the pipe and returns the
            %heat loss rate in kW/m
            %setting up pipe properties.
            ksteel=41; %(W/mK)
            kinsul=0.07; %(W/mK) --> thermal condictivity of calcium silicate 
            Dpipe=0.8128; %pipe outer diameter (m)
            Tpipe=0.015875; %pipe thickness (m)
            Tinsul=0.4064; %insulation thickness (m)
            D=Dpipe+2*Tinsul; %total diameter (m)
            %assume Tf
            Tf=350; %(K)
            kair=0.03; %(W/mK)
            Pr=0.7;
            alpha=29.9*10^-6; %(m^2/s)
            viscosity=20.92*10^-6; %(m^2/s)
            %assume Tinf
            Tinf=303; %(K)
            Beta=1/Tf;
            Ts=2*Tf-Tinf; %(K)
            Ts2=0;
            %creating the iterative loop
            iter=0;
            while abs(Ts-Ts2) > .001;
                Ra=(9.81*Beta*(Ts-Tinf)*D^3)/(viscosity*alpha);
                Nu=(0.6+(0.387*Ra^(1/6))/(1+(0.559/Pr)^(9/16))^(8/27))^2;
                h=(Nu*kair)/D;
                R1=log((Dpipe/2)/((Dpipe-2*Tpipe)/2))/(2*pi()*ksteel); %conductive resistance of pipe (m*K/W)
                R2=log((Dpipe/2+Tinsul)/(Dpipe/2))/(2*pi()*kinsul); %conductive resistance of insul (m*K/W)
                Rconv=1/(h*pi()*D); %(m*K/W)
                Rtotal=R1+R2+Rconv;
                Tdelta=(internal_temp+273)-Tinf; %check the signs of the equations, make sure energy is flowing out t1>tinf
                Q=Tdelta/Rtotal; %W/m
                Ts1=Q*Rconv+Tinf;
                Ts2=Ts;
                Ts=Ts1;
                iter=iter+1;
            end
            Q_loss = Q/1000; %converting from W/m to kW/m
        end  
    end
    
    methods (Access = private)
        function object = setup_arrays(object,discharge)
            if discharge == 1 %discharge on
                object.N = round(object.discharge_time/object.time_step)+1;
            else %discharge off (storage on)
                object.N = round(object.store_time/object.time_step)+1;
            end
            object.time=zeros(object.N,1);
            object.time(1)=0.;
            for count=1:object.N-1
                object.time(count+1)=object.time(count)+object.time_step;
            end
            %initialize variables
            object.p=zeros(object.N,1);  % bar
            object.m1=zeros(object.N,1);  % kg
            object.m2=zeros(object.N,1);
            object.Vol1=zeros(object.N,1);  % m3
            object.Vol2=zeros(object.N,1); %m3
            object.v1=zeros(object.N,1); %m3/kg
            object.v2=zeros(object.N,1); %m3/kg
            object.rho1=zeros(object.N,1);   % kg/m3
            object.rho2=zeros(object.N,1);
            object.h1=zeros(object.N,1);  % kJ/kg 
            object.h2=zeros(object.N,1);
            object.h2in=zeros(object.N,1);
            object.h2out=zeros(object.N,1);
            object.x=zeros(object.N,1);
            object.t1=zeros(object.N,1);   % C
            object.t2=zeros(object.N,1);
            object.m1b=zeros(object.N,1);   % kg/s
            object.m2b=zeros(object.N,1);
            object.mh1b=zeros(object.N,1);  % kJ/s
            object.mh2b=zeros(object.N,1);
            object.mpt1=zeros(object.N,1);   % kg/s
            object.mpt2=zeros(object.N,1);   % kg/s; mpt2=-mpt1 
            object.QLOSS=zeros(object.N,1); % kW, heat loss out of pipe
            object.qloss1=zeros(object.N,1); %kW
            object.qloss2=zeros(object.N,1); %kW
            object.dv1dh=zeros(object.N,1);
            object.dv2dh=zeros(object.N,1);
            object.dv1dp=zeros(object.N,1);
            object.dv2dp=zeros(object.N,1);
            object.Area_int=zeros(object.N,1); %m^2/m^3
            object.Qint=zeros(object.N,1); %kW
            object.VAP_OUT=zeros(object.N,1); %kg/s
            
        end
        function object = initial_conditions(object, X0, P0, VTANK)
         
            object.x(object.i)=X0; %initial quality in tank
            object.p(object.i)=P0; %initial pressure in tank (bar)
            object.rho1(object.i)=XSteam('rhoL_p',object.p(object.i)); %kg/m3
            object.rho2(object.i)=XSteam('rhoV_p',object.p(object.i)); %kg/m3
            object.t1(object.i)=XSteam('Tsat_p',object.p(object.i));
            object.t2(object.i)=object.t1(object.i);
            object.mix_rho=1./(object.x(object.i)/object.rho2(object.i)+(1.-object.x(object.i))/object.rho1(object.i));
            object.m_total=object.mix_rho*VTANK; %total mass in tank
            object.m_original = object.m_total;   % save initial mass inventory
            object.m1(object.i)=object.m_total*(1-object.x(object.i)); %total liquid mass in tank (kg)
            object.m2(object.i)=object.m_total*object.x(object.i); %total vapor mass in tank (kg)
            object.Vol1(object.i)=object.m1(object.i)/object.rho1(object.i); %liquid volume in tank (m3)
            object.Vol2(object.i)=object.m2(object.i)/object.rho2(object.i); %vapor volume in tank (m3)
            object.Vol_total(object.i)=object.Vol1(object.i)+object.Vol2(object.i); %volume of liquid and vapor in tank should be conserved(m3)
            object.v1(object.i)=1/object.rho1(object.i); %liquid specific volume in tank (m3/kg)
            object.v2(object.i)=1/object.rho2(object.i); %vapor specific volume in tank (m3/kg)
            object.h1(object.i) = XSteam('h_pT',object.p(object.i),object.t1(object.i)-object.EPSILON); %(kj/kg)
            object.h2(object.i) = XSteam('h_pT',object.p(object.i),object.t2(object.i)+object.EPSILON);
            object.Dh=(XSteam('v_ph',object.p(object.i),object.h1(object.i)+object.deltaH)-XSteam('v_ph',object.p(object.i),object.h1(object.i)))/abs(object.deltaH); %numerical derivative
            object.Dp=(XSteam('v_ph',object.p(object.i)+object.deltaP,object.h1(object.i))-XSteam('v_ph',object.p(object.i),object.h1(object.i)))/abs(object.deltaP); %numerical derivative
            object.QLOSS(object.i) = object.length*Q_loss(object,(object.t1(object.i)+object.t2(object.i))/2); %initial heat loss value (kW) 
            object.qloss1(object.i)=(object.Vol1(object.i)/(object.Vol1(object.i)+object.Vol2(object.i)))*object.QLOSS(object.i); %volume fraction of heat lost by liquid (kW)
            object.qloss2(object.i)=(object.Vol2(object.i)/(object.Vol1(object.i)+object.Vol2(object.i)))*object.QLOSS(object.i); %volume fraction of heat lost by vapor (kW) 
            object.Area_int(object.i)=2*object.radius*cos(asin(((object.Vol1(object.i)/object.Vol_total(object.i))-object.radius)/object.radius))*object.length; %Minimum Interface Area between the two phases[m^2/m^3]
            object.Qint(object.i)=0; %kW
            object.HH_initial=object.m1(object.i)*object.h1(object.i)+object.m2(object.i)*object.h2(object.i); %kJ
            
        end

        function [md,mu] = rankine(~, po, press)
                %RANKINE A simulation of the ideal Rankine Cycle
        %   RANKINE(IN) generates the thermodynamic properties of the ideal
        %   rankine cycle.  
        %
        %   Calculated parameters:
        %
        %       md    - mass flow rate (kg/s)
        %       mu    - thermodynamic efficiency
        %       bwr   - back work ratio
        %       Qdin  - Rate of energy in (W)
        %       Qdout - Rate of energy out (W)
        %       mdcw  - Condenser mass flow rate
        %
        % This code requires XSteam.m, which available from MATLAB Central.

   
           p_low_side  = 0.07; %bar
           coolant_in = 15; %deg C
           coolant_out = 40; %deg C
         

        %% Look up thermodynamic properties at State 1

           h_location1  = XSteam('hV_p', press)  ;
           s_location1  = XSteam('sV_p', press)  ;

        %% Look up and compute properties for State 2

           s_location2  =  s_location1                  ;

           sg  = XSteam('sV_p', p_low_side)  ;
           sf  = XSteam('sL_p', p_low_side)  ;

           hg  = XSteam('hV_p', p_low_side)  ;
           hf  = XSteam('hL_p', p_low_side)  ;
           hfg = hg - hf              ;

           x_location2  = (s_location2-sf) / (sg-sf)    ;

           h_location2  = hf + x_location2*hfg          ;

        %% State 3 properties

           h_location3  = hf                   ;

        %% State 4 properties

           v_location3  = XSteam('vL_p', p_low_side)  ;
           h_location4  = h_location3 + v_location3*(press - p_low_side)*1e3;

        %% Compute the Thermodynamic Efficiency 

           mu  = ((h_location1-h_location2) - (h_location4-h_location3)) / (h_location1-h_location4);

        %% Compute the Backwork Ratio

           bwr = (h_location4-h_location3) / (h_location1-h_location2) ;

        %% Compute the Mass flow rate at the condenser
           
             md  =  (po*10^6)/( ((h_location1-h_location2)-(h_location4-h_location3))*1000 ) ;
           
        %% Qd (Energy flow)

           Qdin  = md * (h_location1-h_location4)*1000 ; 
           Qdout = md * (h_location2-h_location3)*1000 ;

        %% Steady state energy

           hC    = XSteam('hL_T', coolant_in)  ;
           hH    = XSteam('hL_T', coolant_out)  ;

           mdcw  = md * (h_location2-h_location3)/(hH-hC)      ;

        end    
    end
    
    methods (Static)

        function accumulator = setup_accumulator(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge)
            accumulator = steam_accumulator_separate(T_END,T_STORE,T_RAMP,DT,VTANK,P0,X0,LTANK,RTANK,discharge);
        end
        
        
         
    end
  
end