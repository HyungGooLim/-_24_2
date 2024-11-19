function f_xyz = State_Lorenz(t,f_0)

    global Sigma Rho Beta
    
    dxdt = Sigma*(f_0(2)-f_0(1));
    dydt = f_0(1)*(Rho-f_0(3))-f_0(2);
    dzdt = f_0(1)*f_0(2)-Beta*f_0(3);
    fprintf('time :%d \n',t);
    f_xyz = [dxdt dydt dzdt]';
end