%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-6_1 %%%
%%% 24114529_임형구    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

muE = 3.986004415e5;      % Earth
R = 6378;                 % Earth Radius
r = linspace(R,10*R,100); % Radius of Sphere

% Gravity Potential Function HW-1%
for i = 1:length(r)
    V_0(i) = - muE/r(i);
end

figure();
plot(r,V_0);
ylabel("Gravity potential[km/s^2]");
xlabel("Radius of sphere[km]");
title("Spherical influence case");