close all
clear all
clc

%% Get Kepler Elements %%
global u_e a n
u_e = 3.986004415e5;                     % gravitational coefficient
a = [1e4 1.25e4 0.8696e4];                     % Semi-major axis [km]

ecc = [0 0.2 0.15];
inc = deg2rad(0);                        %[deg]
w =  deg2rad(63.2828);                   %[deg] Argument of perigee
% nu = deg2rad([1:1:361]);               %[rad] True anomaly
simul_time = 27000;                           
RAAN = deg2rad(69.3305);                 %[deg]

%% Kepler 6 Parameters %%
% ecc, w, nu, RAAN, inc, a %

for i=1:3
    n(i) = sqrt(u_e/a(i)^3);                     % mean motion [rad/sec]
    nu(i,:) = n(i)*[1:1:simul_time];             %[rad] True anomaly     
    P = a(i)*(1-ecc(i)^2);
    if i~=3
        for k=1:simul_time
            r(:,k) = P/(1+ecc(i)*cos(nu(i,k)));
            X(i,k) = r(:,k)*cos(nu(i,k));
            Y(i,k) = r(:,k)*sin(nu(i,k));
        end
        if i==1
            plot(X(i,:),Y(i,:),'bo');
            hold on;
        else
            plot(X(i,:),Y(i,:),'ro');
            hold on;
        end
    else
                for k=1:simul_time
            r(:,k) = P/(1+ecc(i)*cos(nu(i,k)));
            X(i,k) = -r(:,k)*cos(nu(i,k));
            Y(i,k) = -r(:,k)*sin(nu(i,k));
        end
        plot(X(i,:),Y(i,:),'co');
        hold on;
    end
end
grid on;
axis equal;
legend('e=0[A]','e=0.2[B]','e=0.15[C]');
ylabel("[km]");
xlabel("[km]");