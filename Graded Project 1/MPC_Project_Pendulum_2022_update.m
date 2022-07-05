% Electrically driven inverted pendulum
% Based on Section 8-5-2 pag 171 from Hanema J., Anticipative MPC for LPV,
% 2018.
% ---------------------------------------------------------------------
% Inthis numerical example is to control the angle q of the electrically
% driven inverted pendulum. Define the state vector
% x = [q dq i] = [ x1 x2 x3] where :
%               q ->  angle of the pendulum   [rad]
%               dq -> angular velocity        [rad/s]  
%               di -> motor current           [A]  
% 
% --------------------------------------------------------------------
clc, clear,
%% ------------------------ Simulation Parameters ---------------------- %%
%           Time
T_step =  1000;     % [s] 
pu = 0.5;           % [fr] step amplitude
T_sample = 0.004;   % [s]
T_sim = 10000;      % []

%% -------------------------- Model Parameters ------------------------- %%

 R = 1;         %   [Ohm] Electrical resistance
 L = 1/1000;    %   [H] Electrical inductance
 k = 6/100;     %   [NA^-1] Motor constant
 b = 1/1000;    %   [Nsm^-1] Friction coefficent
 m = 7/100;     %   [kg] Pendulum mass 
 l = 1/10;      %   [m] Pendulum length
 J = m*l^2;     %   [kgm^2] Pendulum inertia
 g = 9.81;      %   [ms^-2] Standard gravity
 
%% ------------------------ Initial States ------------------------- %%

% For assignment 1. change the initial states in this section.


x1_0 = pi/4;                        % [rad]    q      Angle of the pendulum 

x2_0 = 0;                           % [rad/s]  dq      Angular Velocity 

x3_0 = -(m*g*l*sin(x1_0))/k;      % [A]      i       Motor Current

u_input = x3_0;                     % [V]      u       DC Voltage

%% ---------------------------- Simulation ----------------------------- %%


                    sim('Pendulum_Nonlinear_System')


%% ------------------------------- Plots ------------------------------- %%
figure('Name','System Plot','NumberTitle','off')
set(gca,'FontWeight','bold')
hold on
plot(x1_state,'LineWidth',1);
plot(x2_state,'LineWidth',1);
plot(x3_state,'LineWidth',1);
xlim([0 T_sim]);
xlabel('Time (seconds)');ylabel("States");
legend('$q(t)$','$\dot{q}(t)$','$i(t)$','Interpreter','latex','FontSize',14);
title("System Plot");
hold off


%% ------------------------ Linearised Model --------------------------- %%
% Assignment 2a --> Linearising the Simulink model and convert it to
% discrete system

A_tilde = [0,1,0;
          m*g*l*cos(x1_0)/J, -b/J, k/J;
          0, -k/L, -R/L];
B_tilde = [0;0;R/L];

sysd = c2d(ss(A_tilde,B_tilde,eye(3),0),T_sample);

theta = pi/4;
x_ss = [theta;0;-(m*g*l*sin(theta))/k;];    % Steady States of X
u_ss = x_ss(3);             % Steady State of U

%% ----------------------- Initial Conditions -------------------------- %%
x1_0 = 1.5*pi/4; % changed into 1.5*pi/4 for assignment 2c
x2_0 = 0;   % changed into 0 for assignment 2c
x3_0 = -1.549;  % changed into -1.549 for assignment 2c

x0 = [x1_0;x2_0;x3_0]-x_ss;

%% ----------------- Initial Condition MPC Controller ------------------ %%
% Assignment 2b --> Designing the MPC with constraints
% Assignment 2c --> Q and Re are changed for tuning 

N = 10; % Prediction Horizon
Q = eye(3); Q(1,1) = 65; Q(2,2) = 40; Q(3,3) = 120; % State weight
Re = 160;    % Control input weight

[~,P,~] = dlqr(sysd.A,sysd.B,Q,Re); % Cost function
% P = [31500 1471.70 125; 1471.70 350 30;125 30 125];
% Prediction matrices for MPC
[Phi, Gamma] = ABN2PhiGamma(sysd.A,sysd.B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,Re,P,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

% constrained MPC parameters
xmax = [2*pi;12;8]-x_ss;
xmin = [-(2*pi);-12;-8]-x_ss;
umax = 10-u_ss;
umin = -10-u_ss;
[W, Lu, c] = getWLc(sysd.A, sysd.B, xmax, xmin, umax, umin, Gamma, Phi);

%% -------------------- Closed Loop MPC Controller --------------------- %%
k_sim = 800;
z = zeros(3,k_sim+1);
z(:,1) = x0;
options_qp =  optimoptions('quadprog','Display','off');

x_cloop(:,1) = z(:,1)+x_ss;

for i = 1:k_sim
    warning off
    [u_qp,fval,exitflag] = quadprog(G,F*z(:,i),Lu,c+W*z(:,i),[],[],[],[],[],options_qp);
    if exitflag ~= 1
        warning('exitflag quadprog =%d\n', exitflag)
        if exitflag == -2
            fprintf('Optimization problem is infeasible. \n')
        end
    end
    fprintf("Closed-loop Iteration: %d \n",i);
    u_qp_list(:,i) = u_qp;
    v_mpc(i) = u_qp(1);
    z(:,i+1) = sysd.A*z(:,i) + sysd.B*v_mpc(i); % Closed-loop MPC
    u_mpc_cloop(i) = v_mpc(i) + u_ss;
    x_cloop(:,i+1) = z(:,i+1) + sysd.A*x_ss+ sysd.B*u_ss;
end
%% ---------------------- Simulink MPC Controller ---------------------- %%
clear x1_state x2_state x3_state u_qp u_qp_list exitflag;
T_sim = 1;  % Easier to manage in the loop (speeding the 
k_simulink = 800; % 
x = zeros(3,k_simulink+1);
realX = x;
x(:,1) = x0;
realX(:,1) = x0+x_ss;

options_qp =  optimoptions('quadprog','Display','off');

for i = 1:k_simulink
    [u_qp,fval,exitflag] = quadprog(G,F*x(:,i),Lu,c+W*x(:,i),[],[],[],[],[],options_qp);
    if exitflag ~= 1
        warning('exitflag quadprog =%d\n', exitflag)
        if exitflag == -2
            fprintf('Optimization problem is infeasible. \n')
        end
    end
    fprintf("Simulink Iteration: %d \n",i);
    u_qp_list(:,i) = u_qp;
    u_mpc(i) = u_qp(1); % u_mpc from the quadprog
    realU(i) = u_mpc(i)+u_ss; % Tracking the real u_mpc
    u_input = u_mpc(i)+u_ss; % feeding the u_mpc to the simulink model
    sim('Pendulum_Nonlinear_System');
    x(:,i+1) = [x1_state(2,:)-x_ss(1); x2_state(2,:)-x_ss(2); x3_state(2,:)-x_ss(3)];
    x1_0 = x1_state(2,:); x2_0 = x2_state(2,:); x3_0 = x3_state(2,:); 
    realX(:,i+1) = x(:,i+1)+x_ss;   % Tracking the real states from the Simulink
end

%% ------------------ Closed-loop MPC Controller Plot ------------------ %%
figure('Name','Closed Loop MPC','NumberTitle','off');
subplot(2,1,1);
plot(0:k_sim,x_cloop(1,:),'Linewidth',1);
grid on
hold on
plot(0:k_sim,x_cloop(2,:),'Linewidth',1);
plot(0:k_sim,x_cloop(3,:),'Linewidth',1);
set(gca,'FontWeight','bold')
% xlim([0 k_sim]);ylim([-2 3]);set(gca,'FontWeight','bold')
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',14)
xlabel('$k$','Interpreter','latex','FontSize',14);ylabel('$States$','Interpreter','latex','FontSize',14);
title('State trajectory - Closed-loop Linearized System','FontSize',14);
subplot(2,1,2);
plot(0:k_sim-1,u_mpc_cloop,'Linewidth',1);
grid on
set(gca,'FontWeight','bold')
% % % xlim([0 k_sim]);ylim([-2 0]);set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',14);ylabel('$u$','Interpreter','latex','FontSize',14);
title('Control Input - Closed-loop Linearized System','FontSize',14);

%% ------------------- Simulink MPC Controller Plot -------------------- %%
figure('Name','Simulink Results','NumberTitle','off');
subplot(2,1,1);
plot(0:k_simulink,realX','Linewidth',1);
grid on
set(gca,'FontWeight','bold')
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',14)
xlabel('$k$','Interpreter','latex','FontSize',14);ylabel('$States$','Interpreter','latex','FontSize',14);
title('State trajectory - Non-linear System','FontSize',14);
subplot(2,1,2);
plot(0:k_simulink-1,realU,'Linewidth',1);
grid on
set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',14);ylabel('$u$','Interpreter','latex','FontSize',14);
title('Control Input - Non-linear System','FontSize',14);

