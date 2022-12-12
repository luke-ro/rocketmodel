% Rachael Carraras, Abby Durrell, Kevin Pipich, Luke Roberson, Tyler Schwinck, 
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
% a = 0.047*0.0254;  % [in/s]-> [m/s] burn rate coefficient
a = 6.99e-5; 
n = 0.321; % [-] burn rate exponent, given
cstar = 1500; % [m/s] characteristic velocity
h_grain = 1.505*0.0254; % [in]->[m] motor grain height
r_grain_i = 0.177/2*0.0254; % [in]->[m] motor grain inner radius
r_grain_o = 0.908/2*0.0254; % [in]->[m] motor grain outer radius
r_throat = 0.123/2*0.0254; % [in]->[m] throat radius
r_exit = 0.231/2*0.0254; % [in]->[m] exit radius
Mass = 0.025; % [kg]

%% QUANTITY CALCULATIONS
Vol = h_grain*(r_grain_o^2-r_grain_i^2)*pi(); % [m^3]
rho_p = Mass/Vol; % [kg/m^3]  This is verified
A_throat = pi*(r_throat)^2; % [m^2]
A_exit = pi*(r_exit)^2; % [m^2]
AR_sup = A_exit/A_throat; % supersonic area ratio

%V_burn = 0; % [m^3]
%V_chamber = 0; % [m^3]
j = 1;
% while rb < (r_grain_o - r_grain_i) && rb < h_grain % while there is unburned grain remaining

while true
    
    [A_burn(j), V_burn(j), V_chamber(j)] = burn_geometry(r_grain_i,r_grain_o,h_grain,bd(j)); % [m] burn area, burn cavity volume
    A_burn(j) = A_burn(j) * 0.98;
    Pc(j) = ((a * rho_p * A_burn(j) * cstar) / (A_throat)).^((1)/(1-n)); % [Pa] chamber pressure
    burn_rate(j) = a.*(Pc(j)).^n; % [m/s] burn rate
    bd(j+1) = bd(j) + burn_rate(j) * t_step; % [m] updates burn displacement
    
    % delta_Vol = ; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar,mdot(j),T_metric(j)] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup); %, delta_Vol);
    cstar = cstar*cstar_eff*0.97; % [m/s]
    isp(j) = T_metric(j)./mdot(j)/9.8;
    
    % break condition for the while loop. Makes more sense?
    if A_burn(j) <= 0
        break
    end
   
    
    j = j+1;
end

T_predicted = [0 T_predicted];

t = (0:t_step:t_step*(j-1));

figure
m = 4; n=1;
subplot(m,n,1)
plot(t,A_burn)
title("A\_burn vs t")
xlabel("time")

subplot(m,n,2)
plot(t,V_burn)
title("V\_burn vs t")
xlabel("time")

subplot(m,n,3)
plot(t,V_chamber)
title("V\_chamber vs t")
xlabel("time")

subplot(m,n,4)
plot(t,bd(1:end-1))
title("bd (burn disp) vs t")
xlabel("time")

figure
m = 2; n=1;
subplot(m,n,1)
plot(t,Pc)
title("Chamber pressure (Pc) vs t")
xlabel("time")

subplot(m,n,2)
plot(t,burn_rate)
title("Burn Rate (burn\_rate) vs t")
xlabel("time")

%%
% Experiemental data
%Reading in .csv file
A = readmatrix('data/Static_Fire_24.csv');
T = A(:,1);
t_exp = A(:,2);
t_exp = t_exp-t_exp(1);
t_exp =t_exp./1000;
i = 2;
while i <= length(T)
   if abs(T(i-1)-T(i)) > 0.5
      index = i-2; 
      i = length(T);
   end
   i = i + 1;
end
slope = (T(167)-T(96))/(t_exp(167)-t_exp(96));
T(1:95) = T(1:95)-T(1);
T(166:end) = T(166:end)-T(end);
T(96:167) = T(96:167)-(slope.*t_exp(96:167)-0.6693);

t_exp = t_exp(96:170);
T = T(96:170);
t_exp = t_exp-t_exp(1);
t = (0:t_step:t_step*(j-0));

total_impulse_imp = sum(t_step*T_predicted);
total_impulse_met = sum(t_step*T_metric);

A = readmatrix('Thrust.csv');
t_vendor = A(:,1);
t_vendor = t_vendor-t_vendor(1);
T_vendor =  A(:,2);
t = t(3:end);
t = t-t(1);
T_predicted = T_predicted(3:end);

figure
% plot(t_exp,T,'k','LineWidth',2);
plot(t_vendor,T_vendor,'r','LineWidth',2)
hold on
grid on
grid minor
xlabel('Time [s]')
ylabel('Thrust [lbf]')
plot(t,T_predicted,'b','LineWidth',2)
title("Thrust-Time Profile")
legend('Experimental','Predicted');

% figure
% hold on
% grid on
% xlabel('Time [s]')
% ylabel('ISP [s]')
% plot(t(1:end-3),isp,'LineWidth',2)
% title("ISP vs t")
