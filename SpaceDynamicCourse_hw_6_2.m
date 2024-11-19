%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-6_2 %%%
%%% 24114529_임형구    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% close all;
% clc;

muE = 3.986004415e5;      % Earth
R = 6378;                 % Earth Radius
r = linspace(R,10*R,100); % Radius of Sphere

J2 = 1082.63e-06;
J3 = -2.52e-06;
phi = deg2rad(0:5:90);

figure();
% Gravity Potential Function HW-2%
for j = 1:length(phi)
    for i = 1:length(r)
        zonal_harmonics(j,i) = - (J2/2)*(muE/r(i))*(R/r(i))^2*(3*sin(phi(j))^2-1)-(J3/2)*(muE/r(i))*(R/r(i))^3*(5*sin(phi(j))^3-3*sin(phi(j)));
    end
    plot(r,zonal_harmonics(j,:));
    legends{j} = sprintf("phi: %.1f [deg]", rad2deg(phi(j))); % Store legend label
    hold on;
end

ylabel("Zornal Perturbation Accelerations [km/s^2]");
xlabel("Radius of sphere [km]");
title("Ellipsoidal influence case");
legend(legends); % Set the legend after the loop

figure();
for j = 1:length(phi)
    V_r(j,:) = V_0 + zonal_harmonics(j,:);
    
    plot(r,V_r(j,:));
    legends{j} = sprintf("phi: %.1f [deg]", rad2deg(phi(j))); % Store legend label
    hold on;

end

ylabel("Gravity potential[km/s^2]");
xlabel("Radius of sphere[km]");
title("Ellipsoidal influence case");