function [ThSM_en, Cstar] = thrust_calc(Pa, Pc, Ae, rho_p, burn_rate, A_burn, AR_sup)
    %Pc_en = Pc * 145.038; % [psi] chamber pressure

    %% Aerothermochemistry Data
    ERROR = 0;
    %try % tests if there is any output
        % OUTPUT1 GIVES VALUES IN THE CHAMBER
        % OUTPUT2 GIVES VALUES AT THE NOZZLE THROAT
        % OUTPUT3 GIVES VALUES AT NOZZLE EXIT
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
    
    % Add other aerothermochemistry parameters as needed
        
    %end

    %% MASS FLOW CALCULATION
    m_dot = 0.000005; % [kg/s] propellant mass flow rate

    %% THRUST CALCULATION
    ThSM = 20; % [N] thrust SI
    ThSM_en = ThSM * 0.224809; % [lbf] imperial thrust to match curve data