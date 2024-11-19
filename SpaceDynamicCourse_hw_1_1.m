%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-1 %%%
%%% 24114529_임형구  %%%
%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;


% Lorenz Equation Modeling %
global Sigma Rho Beta
tspan = 0:0.002:0.01;
f_0 = [2 1 1]';
Sigma = 10;
Rho = 28;
Beta = 8/3;
simul_time = 1500;
f = zeros(3,simul_time);

for k = 1:simul_time
    [t,f_xyz] = ode45(@State_Lorenz_hw_1_1,tspan,f_0);
    f(:,k) = f_xyz(end,:);
    f_0 = f(:,k);
    fprintf('step :%d \n',k);
end

% for k = 1:10:simul_time
%     % Update the plot
%     plot3(f(1,k), f(2,k), f(3,k),'-o','Color','b','MarkerSize',2);  % Plot current point
%     hold on;
% 
%     % Setting axis limits
%     xlim([-30 30]);
%     ylim([-30 30]);
%     zlim([0 60]);
% 
%     grid on;
%     title('Lorenz Attractor');
%     % 
%     % % Capture the frame as an image
%     drawnow;
% 
% end

for k = 1:simul_time
    if k == 1
        % 첫 번째 포인트만 있을 때는 이전 점이 없으므로 점만 찍음
        plot3(f(1,k), f(2,k), f(3,k), '-o', 'Color', 'b', 'MarkerSize', 2);  
    else
        % 이전 좌표와 현재 좌표를 연결하여 선을 그리기
        plot3(f(1,k-1:k), f(2,k-1:k), f(3,k-1:k), '-o', 'Color', 'b', 'MarkerSize', 2);
    end
    hold on;
    
    % Setting axis limits
    xlim([-30 30]);
    ylim([-30 30]);
    zlim([0 60]);
    
    grid on;
    title('Lorenz Attractor');
    
    % Capture the frame as an image
    drawnow;
end