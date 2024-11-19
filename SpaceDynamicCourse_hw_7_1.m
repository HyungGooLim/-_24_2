%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-7_1  %%%
%%% 24114529_임형구     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all;
clc;

% 물리 상수
global mu
mu = 398600;  % km^3/s^2 (지구 중력 상수)
J2 = 0.00108263;
R = 6378;                 % Earth Radius

% 궤도 요소
a = 7001;                  % 반장축 [km]
ecc = 0.01;                  % 이심률
inc = 45 * pi / 180;  % 궤도 경사 [rad]
Omega = 45 * pi / 180 + (-1.5*J2*n*());      % 상승노드 경도 [rad]
w = 45 * pi / 180;          % 근점 인수 [rad]
M0 = 45 * pi / 180;         % 평균 근점이각 [rad]
h = sqrt(mu*a*(1-ecc^2));   % 각운동량

% 평균운동 계산
n = sqrt(mu/a^3);         % mean motion
tspan = 2*pi*sqrt(a^3/mu);


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
function dxdt = orbitEquations(~, x)
    global mu
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);  % 거리 계산
    dxdt = zeros(6,1);  % 미분 방정식 배열

    % 위치 변화율 (속도)
    dxdt(1) = x(4);
    dxdt(2) = x(5);
    dxdt(3) = x(6);

    % 속도 변화율 (가속도)
    dxdt(4) = -mu * x(1) / r^3;
    dxdt(5) = -mu * x(2) / r^3;
    dxdt(6) = -mu * x(3) / r^3;
end