function [Phi, Gamma] = LPVPhiGamma(B,N,P)
%ABN2PhiGamma Summary of this function goes here
%   Detailed explanation goes here

nx = size(B,1);
nu = size(B,2);

R = 1;         %   [Ohm] Electrical resistance
L = 1/1000;    %   [H] Electrical inductance
k_param = 6/100;     %   [NA^-1] Motor constant
b = 1/1000;    %   [Nsm^-1] Friction coefficent
m = 7/100;     %   [kg] Pendulum mass 
l = 1/10;      %   [m] Pendulum length
J = m*l^2;     %   [kgm^2] Pendulum inertia
g = 9.81;      %   [ms^-2] Standard gravity
Ts = 0.004;

% sdpvar A;
% A{num} = [1,Ts,0; Ts*m*g*l*P(num)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)];

% Phi
Phi = [1,Ts,0; Ts*m*g*l*P(1)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)];
if N>1
        for j = 2:N
            temp = [1,Ts,0; Ts*m*g*l*P(j)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)];
            for k=j-1:-1:1
                temp = temp*[1,Ts,0; Ts*m*g*l*P(k)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)];
            end
            Phi = [Phi;temp];
        end
        
end

temp = [];
% Gamma
Gamma = zeros(nx*N,nu*N);

for i = 1:N
    % A = [1,Ts,0; Ts*m*g*l*P(i)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)];
    for j = 1:i
        if j == i
            Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = B;
        elseif i>j
            temp = [1,Ts,0; Ts*m*g*l*P(i-1)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)];
            for k=j:1:i-2
                temp = [1,Ts,0; Ts*m*g*l*P(k)/J, 1-(b*Ts/J), k_param*Ts/J; 0 -k_param*Ts/L 1-(Ts*R/L)]*temp;
            end
            Gamma((i-1)*nx+1:nx*i,(j-1)*nu+1:j*nu) = temp*B;
        end
    end
end
end

