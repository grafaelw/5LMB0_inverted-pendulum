clear, clc,
%% -------------------- Simulation & Model Parameters ------------------ %%
%           Time
T_step =  1000;     % [s] 
pu = 0.5;           % [fr] step amplitude
T_sample = 0.004;   % [s]
T_sim = 1;          % []
k_sim = 3000;       % [] Quadprog looping simulation

R = 1;         %   [Ohm] Electrical resistance
L = 1/1000;    %   [H] Electrical inductance
k = 6/100;     %   [NA^-1] Motor constanttem
b = 1/1000;    %   [Nsm^-1] Friction coefficent
m = 7/100;     %   [kg] Pendulum mass 
l = 1/10;      %   [m] Pendulum length
J = m*l^2;     %   [kgm^2] Pendulum inertia
g = 9.81;      %   [ms^-2] Standard gravity

%% ---------------- Linearised Model & Initial States ------------------ %%
% Assignment 1a --> Linearising the Simulink model and convert it to discrete system

A_tilde = [0,1,0;
          m*g*l*cos(pi/4)/J, -b/J, k/J;
          0, -k/L, -R/L];
B_tilde = [0;0;R/L];

sysd = c2d(ss(A_tilde,B_tilde,eye(3),0),T_sample);

theta = pi/4;
x_ss = [theta;0;-(m*g*l*sin(theta))/k;];    % Steady States of X
u_ss = x_ss(3);                             % Steady State of U

%% ----------------------- Initial Conditions -------------------------- %%
x1_0 = 1.1*pi/2; 
x2_0 = 10;   
x3_0 = -1.549;  

x0 = [x1_0;x2_0;x3_0]-x_ss;

%% ---------------- Inactive Constraint MPC Controller ----------------- %%

% Inactive Constraints
N = 10; % Prediction Horizon
Q = eye(3); Q(1,1) = 120; Q(2,2) = 10; Q(3,3) = 10; % State weight
Re = 10;    % Control input weight
P = Q;

[K,~,~] = dlqr(sysd.A,sysd.B,Q,Re); % Terminal cost function & LQR gain

%% ----------------- Initial Condition MPC Controller ------------------ %%

% Prediction matrices for MPC
[Phi, Gamma] = ABN2PhiGamma(sysd.A,sysd.B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,Re,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

% Constrained MPC parameters
xmax = [2*pi;12;8]-x_ss;
xmin = [-2*pi;-12;-8]-x_ss;
umax = 10-u_ss;
umin = -10-u_ss;

[W, Lu, c] = getWLc(sysd.B, xmax, xmin, umax, umin, Gamma, Phi);


%% -------------------- Closed Loop MPC Controller --------------------- %%

z = zeros(3,2);
predX = z;
z(:,1) = x0;
% predX(:,1) = x0+x_ss;
options_qp =  optimoptions('quadprog','Display','off');


warning off
%calculate rho >>> 

[u_qp,fval,exitflag] = quadprog(G,F*z(:,1),Lu,c+W*z(:,1),[],[],[],[],[],options_qp);
if exitflag ~= 1
    warning('exitflag quadprog =%d\n', exitflag)
    if exitflag == -2
        fprintf('Optimization problem is infeasible. \n')
    end
end
% calculating u and x
z(:,2) = sysd.A*z(:,1) + sysd.B*u_qp(1); % Closed-loop MPC
predX(:,1) = z(:,2)+sysd.A*x_ss+sysd.B*u_ss;
for j = 2:N
    predX(:,j) = sysd.A*z(:,1)+sysd.B*u_qp(j);
    predX(:,j) = predX(:,j)+sysd.A*x_ss+sysd.B*u_ss;
end

Rho = sinc(predX(1,1)/pi);
for i=2:N
    Rho = [Rho sinc(predX(1,i)/pi)];
end

Uk = u_qp;

save("Uk-en-Rho.mat","Uk","Rho")