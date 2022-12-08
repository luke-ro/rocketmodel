%% Rocket Project, Analytical Portion
%Housekeeping
clc
clear
close all
%% Problem 1
%Reading in .csv file
A = readmatrix('data/Thrust.csv');
t = A(1:80,1);
T = A(1:80,2); %erroneous data past index 80
plot(t, T)
title("Thrust vs time for a F26T curve")
xlabel('Time [s]')
ylabel('Thrust [lbf]')

% Verify total impulse from predicted thrust curve is close to claimed value in data
%Insert Thrust-impulse relation equation here: Impulse = F*time
F = mean(A(:,2));
deltat = A(end,1);
Impulse = F*deltat*4.44822162%Converting from lb-s to N-s
disp("Claimed Impulse = 50 N-s")
total_impulse_dig = sum((t(2:end)-t(1:end-1)).*T(2:end));

%% remove bias from thrust Thrust
%Reading in .csv file
A = readmatrix('data/Static Fire 24.csv');
T = A(:,1);
t = A(:,2);
t = t-t(1);
t  =t./1000;
i = 2;
while i <= length(T)
   if abs(T(i-1)-T(i)) > 0.5
      index = i-2; 
      i = length(T);
   end
   i = i + 1;
end
slope = (T(167)-T(96))/(t(167)-t(96));
T(1:95) = T(1:95)-T(1);
T(166:end) = T(166:end)-T(end);
T(96:167) = T(96:167)-(slope.*t(96:167)-0.6693);



plot(t,T);
grid on
title("Thrust vs time for static fire 24")
xlabel('Time [s]')
ylabel('Thrust [lbf]')

%% Problem 2: Plugging in Values to get Ab-- done in burn_geometry

%% Problem 3: Develop aerothermochemistry tables or algorithms for determination of these parameters.  
% Provide evidence supporting your choice of ingredients and percent mass for the propellant
% composition. 


