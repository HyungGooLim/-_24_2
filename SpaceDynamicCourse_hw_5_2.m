%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-5_2 %%%
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
xf = [-1.764e02, 6.206e03, -2,735e03, -5.250, 2.124e0, 5.167]';

% ODE 설정 옵션: 상대/절대 오차 설정
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% 시간 배열 생성
tspan = linspace(t0, tf, 1000);  % 1000개 시간 지점

%% Solver 1: ode45 실행 및 시간 측정
tic;  % 시간 측정 시작
[t_45, sol_45] = ode45(@orbitEquations_J2, tspan, x0, options);
time_45 = toc;  % 시간 측정 종료

%% Solver 2: ode113 실행 및 시간 측정
tic;  % 시간 측정 시작
[t_113, sol_113] = ode113(@orbitEquations_J2, tspan, x0, options);
time_113 = toc;  % 시간 측정 종료

%% 결과 출력 및 비교
fprintf('ode45 계산 시간: %.6f 초\n', time_45);
fprintf('ode113 계산 시간: %.6f 초\n', time_113);

% 결과 시각화 ode45 (3D 궤도 플롯)
figure;
plot3(sol_45(:,1), sol_45(:,2), sol_45(:,3), 'r', 'LineWidth', 1);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit Trajectory - ODE45');
grid on;
axis equal;

% 결과 시각화 ode113 (3D 궤도 플롯)
figure;
plot3(sol_113(:,1), sol_113(:,2), sol_113(:,3), 'r', 'LineWidth', 1);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit Trajectory - ODE113');
grid on;
axis equal;
% 
% hold on;
% 
% load('topo.mat','topo','topomap1');
% Re=6378.137;
% 
% phi=linspace(0,pi,30);
% theta=linspace(0,2*pi,40);
% [phi,theta]=meshgrid(phi,theta);
% xsphere=Re*sin(phi).*cos(theta);
% ysphere=Re*sin(phi).*sin(theta);
% zsphere=Re*cos(phi);
% mhndl1=mesh(xsphere,ysphere,zsphere);
% h =surface(xsphere,ysphere,zsphere,'FaceColor','texture','CData',topo);
% axis equal
