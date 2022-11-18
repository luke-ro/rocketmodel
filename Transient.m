% Kevin Pipich, Luke Roberson
% ASEN 4013 Rocket Project
close all;clear;clc;

% SI UNITS used for this code.  Conversions from Imperial to SI system made
% as appropriate

t(1) = 0; % [s] initial time
bd(1) = 0; % [m] initial burn grain displacement
T_vendor = csvread('Thrust.csv');  % vendor provided data
T_exp = csvread('Test.csv');  % data from CU experimental static fire

%% INPUTS
cstar_eff = 1.; % [-], cstar efficiency
t_step = 0.05; % [s] time step
P_atm = 101325.; % [Pa] ambient pressure
a = .000005; % [-] burn rate coefficient
n = 0.5; % [-] burn rate exponent
cstar = 1500.; % [m/s] characteristic velocity
h_grain = 1.505*0.0254; % [in]->[m] motor grain height
r_grain_i = 0.177/2*0.0254; % [in]->[m] motor grain inner radius
r_grain_o = 0.908/2*0.0254; % [in]->[m] motor grain outer radius
r_throat = 0.123/2*0.0254; % [in]->[m] throat radius
r_exit = 0.231/2*0.0254; % [in]->[m] exit radius
Mass = 0.025; % [kg]

%% QUANTITY CALCULATIONS
Vol = h_grain*(r_grain_o^2-r_grain_i^2)*pi(); % [m^3]
rho_p = Mass/Vol; % [kg/m^3]
A_throat = pi*(r_throat)^2; % [m^2]
A_exit = pi*(r_exit)^2; % [m^2]
AR_sup = A_exit/A_throat; % supersonic area ratio

%V_burn = 0; % [m^3]
%V_chamber = 0; % [m^3]
j = 1;
% while rb < (r_grain_o - r_grain_i) && rb < h_grain % while there is unburned grain remaining

while true
    
    [A_burn(j), V_burn(j), V_chamber(j)] = burn_geometry(r_grain_i,r_grain_o,h_grain,bd(j)); % [m] burn area, burn cavity volume
    
    Pc(j) = ((a * rho_p * A_burn(j) * cstar) / (A_throat)).^((1)/(1-n)); % [Pa] chamber pressure
    burn_rate(j) = a*(Pc(j))^n; % [m/s] burn rate
    bd(j+1) = bd(j) + burn_rate(j) * t_step; % [m] updates burn displacement
    
    % delta_Vol = ; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup); %, delta_Vol);
    cstar = cstar*cstar_eff; % [m/s]
    
    % break condition for the while loop. Makes more sense?
    if A_burn(j) <= 0
        break
    end
    
    j = j+1;
end

t = (0:t_step:t_step*(j-1));

figure
m = 4; n=1;
subplot(m,n,1)
plot(t,A_burn)
title("A\_burn vs t")

subplot(m,n,2)
plot(t,V_burn)
title("V\_burn vs t")

subplot(m,n,3)
plot(t,V_chamber)
title("V\_chamber vs t")

subplot(m,n,4)
plot(t,bd(1:end-1))
title("bd (burn disp) vs t")

figure
m = 3; n=1;
subplot(m,n,1)
plot(t,T_predicted)
title("Thrust vs t")

subplot(m,n,2)
plot(t,Pc)
title("Chamber pressure (Pc) vs t")

subplot(m,n,3)
plot(t,burn_rate)
title("Burn Rate (burn_rate) vs t")









