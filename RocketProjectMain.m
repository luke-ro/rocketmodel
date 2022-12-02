%% Rocket Project, Analytical Portion
%Housekeeping
clc
clear
close all
%% Problem 1
%Reading in .csv file
A = readmatrix('Thrust.csv');
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

%% Problem 2: Plugging in Values to get Ab-- done in burn_geometry

%% Problem 3: Develop aerothermochemistry tables or algorithms for determination of these parameters.  
% Provide evidence supporting your choice of ingredients and percent mass for the propellant
% composition. 


