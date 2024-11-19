%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-7_2  %%%
%%% 24114529_임형구     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
% close all;
clc;

% 물리 상수
global mu J2 R
mu = 398600;  % km^3/s^2 (지구 중력 상수)
J2 = 0.00108263;
R = 6378;  % Earth Radius

% 궤도 요소
a = 7001;  % 반장축 [km]
ecc = 0.01;  % 이심률
P = a * (1 - ecc^2);
inc = 45 * pi / 180;  % 궤도 경사 [rad]
n = sqrt(mu / a^3);  % mean motion
tspan = 2 * pi * sqrt(a^3 / mu);
Omega = 45 * pi / 180;  % 상승노드 경도 [rad]
w = 45 * pi / 180;  % 근점 인수 [rad]
M0 = 45 * pi / 180;  % 평균 근점이각 [rad]
h = sqrt(mu * a * (1 - ecc^2));  % 각운동량
   

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
% ODE 시스템 정의 (운동 방정식)
function dxdt = orbitEquations(t, x)
    global mu J2 R
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2); % 거리 계산
    dxdt = zeros(6,1);  % 미분 방정식 배열
    % 위치에 대 한 변화율
    dxdt(1) = x(4);  % vx
    dxdt(2) = x(5);  % vy
    dxdt(3) = x(6);  % vz
    % 속도에 대한 변화율 (중력 효과)
    dxdt(4) = (-mu*x(1)/r^3)*(1+1.5*J2*(R/r)^2*(1-5*(x(3)/r)^2));  % ax
    dxdt(5) = (-mu*x(2)/r^3)*(1+1.5*J2*(R/r)^2*(1-5*(x(3)/r)^2));  % ay
    dxdt(6) = (-mu*x(3)/r^3)*(1+1.5*J2*(R/r)^2*(3-5*(x(3)/r)^2));  % az
end


% Mean anomaly to True anomaly
M = M0;
E = Mean2Eccen(M, ecc);
nu = atan2(sqrt(1 - ecc^2) * sin(E), cos(E) - ecc);

% Get Distance Satellite from Elliptical Focus
p = a * (1 - ecc^2);
r_0 = p / (1 + ecc * cos(nu));

% Coordinate Transformation
r_ECI = r_0 .* R_pqw_to_eci(Omega, w + nu, inc);
v_ECI = (-mu / h) .* V_pqw_to_eci(Omega, w + nu, w, inc, ecc);

% Get Initial Condition [x0, y0, z0, vx0, vy0, vz0]
x0 = [r_ECI; v_ECI];
x0 = x0(:);  % 벡터 형태로 변환

% Numerical Integration (ODE45)
[t, sol] = ode45(@orbitEquations, [0, tspan*10], x0);

% 결과 시각화 (3D 궤도 플롯)
figure();
plot3(sol(:, 1), sol(:, 2), sol(:, 3), 'r-', 'LineWidth', 0.7);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Orbit Trajectory');
grid on;
axis equal;
hold on;

% 지구 그리기
load('topo.mat', 'topo', 'topomap1');
Re = 6378;

phi = linspace(0, pi, 30);
theta = linspace(0, 2 * pi, 40);
[phi, theta] = meshgrid(phi, theta);
xsphere = Re * sin(phi) .* cos(theta);
ysphere = Re * sin(phi) .* sin(theta);
zsphere = Re * cos(phi);
mhndl1 = mesh(xsphere, ysphere, zsphere);
h = surface(xsphere, ysphere, zsphere, 'FaceColor', 'texture', 'CData', topo);
axis equal
