%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 우주역학특론 HW-8_1  %%%
%%% 24114529_임형구     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;

% 초기조건
% xi= [1.0;0; -1.0;0; 0;sqrt(3); 0;0; 0;0; 0;0;];

% xi= [1;0; -1;0; 0;1; 0;0; 0;0; 0;0]; %case2
 % xi= [1;0.1; -1;0; 0;sqrt(3); 0;0; 0;0; 0;0]; %case3
 xi= [1;0; -1;0; 0;0; 0;1.0; 0;1.0; 0;0]; %case4

global m1 m2 m3
m1 = 1; m2 = 1; m3 = 100;

% Three-Body Problem
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,Y]=ode45(@Three_Body,[1 10],xi,options);
% [t,Y]=ode113(@TBP,[1 10],xi,options);

figure;
plot(Y(:,1),Y(:,2),Linewidth=1.5); hold on;
plot(Y(:,3),Y(:,4),Linewidth=1.5);
plot(Y(:,5),Y(:,6),Linewidth=1.5);
xlabel('x'); ylabel('y')
grid on