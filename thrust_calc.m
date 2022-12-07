function [ThSM_en, Cstar, m_dot] = thrust_calc(Pa, Pc, Ae, rho_p, burn_rate, A_burn, AR_sup)
    % Pa: atmospheric pressure
    % Pc: chamber pressure
    % Ae: area exit
    % rho_p: density of fuel after burn
    % burn_rate: burn rate
    % A_burn: area of the burn
    % AR_sup: Supersonic area ratio 
    
    %Pc_en = Pc * 145.038; % [psi] chamber pressure
    
    % mdot_e = mdot_p + d/dt(rho*V); % mass flow depending on changing
    % density and volume. from cons of mass
    
    % the following
    pc_lookup = [1:1:10].*1000000; %pa
%     pe_lookup = [0.14706 0.29412 0.44118 0.58824 0.73529 0.88235 1.0294 1.1765 1.3235 1.4706].*100000; %Pa
%     ve_lookup = [819.4*3.226 818.3*3.234 817.7*3.238 817.4*3.240 817.1*3.242 817.0*3.243 816.8*3.244 816.7*3.245 816.6*3.246 816.5*3.246];
    pe_lookup = [0.47819 0.94684 1.4133 1.8786 2.3432 2.8073 3.2709 3.7342 4.1972 4.6600 ].*100000;
    ve_lookup = [924.0*2.554 922.4*2.566 921.5*2.572 921.0*2.575 920.6*2.578 920.3*2.580 920.1*2.582 919.9*2.583 919.7*2.584 919.6*2.585];
    %% Aerothermochemistry Data
    ERROR = 0;
    %try % tests if there is any output
        % OUTPUT1 GIVES VALUES IN THE CHAMBER
        % OUTPUT2 GIVES VALUES AT THE NOZZLE THROAT
        % OUTPUT3 GIVES VALUES AT NOZZLE EXIT
        A_throat = Ae/AR_sup;
        [Output] = aerothermochemistry(Pc, AR_sup);
    %catch
    %   ERROR = 1;
    %end
    
    %if ERROR == 1 % sets Mach and alpha to zero if output DNE
    %    alpha = 0;
    %    Mach = 0;
    %    Pe = Pa;
    %    Cstar = 0;
    %    rho_g = 0;
    %else
        Tc = Output.T; % [K] chamber temperature
        gamma_c = Output.g; % ratio of specific heats
        MolWt_c = Output.MolWt; % molecular weight
        rho_c = Output.rho; % [kg/m3] gas density
        Cstar = Output.c; % [m/s] c* velocity
        
        % From CEA, probably wrong
        V_e = interp1(pc_lookup,ve_lookup,Pc,'linear','extrap');
        Pe =  interp1(pc_lookup,pe_lookup,Pc,'linear','extrap');
    % Add other aerothermochemistry parameters as needed   
    
    %end

    %% MASS FLOW CALCULATION
    % m_dot = 0.000005; % [kg/s] propellant mass flow rate
    m_dot = (rho_p-rho_c)*burn_rate*A_burn;         % [kg/s] mass flow
    Cstar = (A_throat*Pc)/m_dot;
    %% THRUST CALCULATION   
    % get V_e and P_e from CEA somehow
    ThSM = ((m_dot*V_e))+((Pe-Pa)*Ae);         % [N] thrust SI
    % ThSM = 20; % [N] thrust SI
    ThSM_en = ThSM * 0.224809; % [lbf] imperial thrust to match curve data
end