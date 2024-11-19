%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-5_1 %%%
%%% 24114529_임형구    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;


global muE J2 R

muE = 3.986004415e5; % Earth
R = 6378;            % Earth Radius
G = 6.67259e-20;     % Gravitational Constant
hours = 3600;        % Conversion variable between seconds & hours
days = hours*24;     % Conversion variable between seconds % days
J2 = 0.00108263;

span = 10;           % days to propagate
t0 = 0;
tf = span*days;      % initail and final times
x0 = [4.803e03, 1.228e03, -5.097e03, -3.87, 6.37e0, -1.73]';    %retrieve initial ephemeris as inital states

% ODE45 solver 호출
tspan = linspace(t0, tf, 10000);  % 1000개의 시간 지점 생성
[t, sol] = ode45(@orbitEquations_J2, tspan, x0);


% 결과 시각화 ode45 (3D 궤도 플롯)
figure;
plot3(sol(:,1), sol(:,2), sol(:,3), 'r', 'LineWidth', 1);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit Trajectory - ODE45');
grid on;
axis equal;

hold on;
load('topo.mat','topo','topomap1');
Re=6378;

phi=linspace(0,pi,30);
theta=linspace(0,2*pi,40);
[phi,theta]=meshgrid(phi,theta);
xsphere=Re*sin(phi).*cos(theta);
ysphere=Re*sin(phi).*sin(theta);
zsphere=Re*cos(phi);
mhndl1=mesh(xsphere,ysphere,zsphere);
h =surface(xsphere,ysphere,zsphere,'FaceColor','texture','CData',topo);
axis equal

grid on;
title('J2 Perturbation');

% 결과 시각화 ode45 (3D 궤도 플롯)
figure;
plot3(sol(:,1), sol(:,2), sol(:,3), 'r', 'LineWidth', 1);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit Trajectory - ODE45');
grid on;
axis equal;