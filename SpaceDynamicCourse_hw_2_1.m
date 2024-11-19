close all;
clear all;
clc;

% ODE 시스템 정의 (운동 방정식)
function dxdt = orbitEquations(t, x)
    global mu
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2); % 거리 계산
    dxdt = zeros(6,1);  % 미분 방정식 배열
    % 위치에 대 한 변화율
    dxdt(1) = x(4);  % vx
    dxdt(2) = x(5);  % vy
    dxdt(3) = x(6);  % vz
    % 속도에 대한 변화율 (중력 효과)
    dxdt(4) = -mu * x(1) / r^3;  % ax
    dxdt(5) = -mu * x(2) / r^3;  % ay
    dxdt(6) = -mu * x(3) / r^3;  % az

    % dxdt = zeros(4,1);  % 미분 방정식 배열
    % % % 위치에 대 한 변화율
    % dxdt(1) = x(3);  % vx
    % dxdt(2) = x(4);  % vy
    % % % 속도에 대한 변화율 (중력 효과)
    % dxdt(3) = x(1)*x(4)^2-x(1)^(-2);
    % dxdt(4) = -x(3)*x(4)/x(1);
end


%Initial Condition%
for k=1:3
    if k==1
        x0 = [9000 0 6000 0.0 5.0 5.0];  % [x, y, z, vx, vy, vz]
        dt = [0:100:86400];  % 시간 간격
    elseif k==2
        x0=[7000 9000 0 -5.0 7.0 0.0];
        dt=[0:100:30000];
    else
        x0=[9000 0 6000 0.0 5.0 5.0]; % ad_x = +1e-5, +1e-6
        dt=[0:100:5*86400];
    end
 

    global mu
    % 중력 상수 및 중심 천체 질량 (예: 지구 기준)
    mu = 398600;  % km^3/s^2 (지구 기준 중력 상수)
    
    
    % ODE45 solver 호출
    [t, sol] = ode45(@orbitEquations, dt, x0);
    
    % 결과 시각화 (3D 궤도 플롯)
    figure;
    plot3(sol(:,1), sol(:,2), sol(:,3), 'LineWidth', 2);
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('3D Orbit Trajectory');
    grid on;
    axis equal;

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
end
