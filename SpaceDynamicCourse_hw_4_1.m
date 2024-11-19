%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-4-2 %%%
%%% 24114529_임형구    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Classical Two Body Problem (CTBP) %

clear all;
close all;
clc;
d2r = pi/180;
mu=398600.5;                                                                %Earth gravitational constant km^3/s^2

% Orbital Elements
a = 7000;                 % Semi major axis
ecc=0.01;                   % eccentricity
inclination=45*d2r;                 % inclination [deg]
Omega = 45*d2r;               % Longitude of the ascending node [deg]
w = 45*d2r;               % argument of perigee [deg]
M0 = 45*d2r;              % mean anomaly at apogee [deg]

n = sqrt(mu/a^3);         % mean motion
h = sqrt(mu*a*(1-ecc^2));   % 각운동량
tspan = 2*pi*sqrt(a^3/mu);
t=0:0.01:tspan*5;


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

%Transformation matrix Definition
function Rotation = C_1(phi)
    C_1 = [1 0 0;...
           0 cos(phi) sin(phi);...
           0 -sin(phi) cos(phi)];
    Rotation = C_1;
end

function Rotation = C_3(psi)
    C_3 = [cos(psi) sin(psi) 0;...
           -sin(psi) cos(psi) 0;...
           0 0 1];
    Rotation = C_3;
end

%Rotation matrix Definition
function Rotation = R_pqw_to_eci(Omega,theta,inc)
    R = [cos(Omega)*cos(theta)-sin(Omega)*sin(theta)*cos(inc);...
         sin(Omega)*cos(theta)+cos(Omega)*sin(theta)*cos(inc);...
         sin(theta)*sin(inc)];
    Rotation = R;
end

function Rotation = V_pqw_to_eci(Omega,theta,w,inc,ecc)
    R = [cos(Omega)*(sin(theta)+ecc*sin(w))+sin(Omega)*(cos(theta)+ecc*cos(w))*cos(inc);...
         sin(Omega)*(sin(theta)+ecc*sin(w))-cos(Omega)*(cos(theta)+ecc*cos(w))*cos(inc);...
         -(cos(theta)+ecc*cos(w))*sin(inc)];
    Rotation = R;
end


for i=1:1:length(t)
    if i==1
        M(1) = M0;
    else
        M(i) = M(i-1)+n*(t(i)-t(i-1));
    end
    E(i) = Mean2Eccen(M(i),ecc);
    nu(i) = atan2((sqrt(1-ecc^2)*sin(E(i))/1-ecc*cos(E(i))),(cos(E(i))-ecc)/(1-ecc*cos(E(i))));
    p(i) = a*(1-ecc^2);
    r_0(i) = p(i) / (1 + ecc*cos(nu(i)));
    X(i) = r_0(i)*cos(nu(i));
    Y(i) = r_0(i)*sin(nu(i));
    r_PQF(:,i) = [X(i) Y(i) 0]';
    % Coordinate Transformation
    r_ECI(:,i) = C_3(Omega)'*C_1(inclination)'*C_3(w)'*r_PQF(:,i);
    % r_ECI(:,i) = r_0(i).*R_pqw_to_eci(Omega,nu(i),inclination);
    % v_ECI(:,i) = (-mu/h).*V_pqw_to_eci(Omega,nu(i),w,inclination,ecc);
end

plot3(r_ECI(1,1:1:length(t)),r_ECI(2,1:1:length(t)),r_ECI(3,1:1:length(t)),'r-', "LineWidth",0.7);
hold on;

load('topo.mat','topo','topomap1');
Re=6378.137;

phi=linspace(0,pi,30);
theta=linspace(0,2*pi,40);
[phi,theta]=meshgrid(phi,theta);
xsphere=Re*sin(phi).*cos(theta);
ysphere=Re*sin(phi).*sin(theta);
zsphere=Re*cos(phi);
mhndl1=mesh(xsphere,ysphere,zsphere);
h =surface(xsphere,ysphere,zsphere,'FaceColor','texture','CData',topo);
axis equal

hold on;
grid on;