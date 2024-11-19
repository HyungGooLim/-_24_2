%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 우주역학특론 HW3 %%%%
%%%% 24114529 임형구 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;


% Parameter %
a = 22500;               % [km]
e = [0:0.2:1];           % eccentricity
mu = 3.986004418e5;      % [km^3/s^-2]earth gravitational constant
n = sqrt(mu/a^3);        % [rad/s]mean motion
M = linspace(-pi,pi,20);


function E = Mean2Eccen(M, e)
    E_n1 = M;
    f = (M-E_n1+e*sin(E_n1));
    fdot= (e*cos(E_n1)-1);
    E_n2 = E_n1-(f/fdot);
    % E_up = M+e*sin(0.1);
    % E_down = E_up + 0.1;
    while (abs(E_n1 - E_n2) > 0.001)
        E_n1 = E_n2;
        f = (M-E_n1+e*sin(E_n1));
        fdot= (e*cos(E_n1)-1);
        E_n2 = E_n1-(f/fdot);
    end 
    E = E_n2;
end


for i=1:1:length(e)
    for j=1:1:length(M)
        E(i,j) = Mean2Eccen(M(j),e(i));
    end
    plot(E(i,:),M);
    % legend(sprintf('ecc = %d', e(i)));
    hold on;
end

hold on;
ylabel("Eccentricity anomaly(E) [rad]");
xlabel("Mean anomaly(M) [rad]");