function dxdt = orbitEquations_J2(t, x)
    global muE J2 R
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2); % 거리 계산
    dxdt = zeros(6,1);  % 미분 방정식 배열
    % 위치에 대 한 변화율
    dxdt(1) = x(4);  % vx
    dxdt(2) = x(5);  % vy
    dxdt(3) = x(6);  % vz
    % 속도에 대한 변화율 (중력 효과)
    dxdt(4) = (-muE*x(1)/r^3)*(1+1.5*J2*(R/r)^2*(1-5*(x(3)/r)^2));  % ax
    dxdt(5) = (-muE*x(2)/r^3)*(1+1.5*J2*(R/r)^2*(1-5*(x(3)/r)^2));  % ay
    dxdt(6) = (-muE*x(3)/r^3)*(1+1.5*J2*(R/r)^2*(3-5*(x(3)/r)^2));  % az

end