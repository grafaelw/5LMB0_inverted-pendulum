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
k = 6/100;     %   [NA^-1] Motor constant
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
x1_0 = 1.1*pi/4; 
x2_0 = 6;   
x3_0 = -1;  

x0 = [x1_0;x2_0;x3_0]-x_ss;

%% ---------------- Inactive Constraint MPC Controller ----------------- %%

% Inactive Constraints
N = 10; % Prediction Horizon
Q = eye(3); Q(1,1) = 120; Q(2,2) = 10; Q(3,3) = 10; % State weight
Re = 10;    % Control input weight

[K,P,~] = dlqr(sysd.A,sysd.B,Q,Re); % Terminal cost function & LQR gain
K_lqr = -K;
% A_cl = sysd.A+sysd.B*K_lqr;

%% -------------- LMI Terminal Cost & Stability Guarantee -------------- %%
O = sdpvar(3,3);Y = sdpvar(1,3);
Z = Q+K_lqr'*R*K_lqr;
% L1 = [O, (sysd.A*O+sysd.B*Y)', O, Y';
%       (sysd.A*O+sysd.B*Y), O, zeros(size(O)), zeros(size(Y'));
%       O, zeros(size(O)), inv(Q), zeros(size(Y'));
%       Y, zeros(size(Y)), zeros(size(Y)), inv(R)]>=0;
L1 = [O,(sysd.A*O+sysd.B*Y)',O;           % Lyapunov's Quadratic Formula
      (sysd.A*O+sysd.B*Y),O,zeros(size(O));
      O,zeros(size(O)),Z^-1]>=0;
L2 = O>=1e-9;
constraints = L1 + L2;
diagnostics = optimize(constraints);
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
P = value(O)^-1;
K = value(Y)*value(O)^-1;

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

%% ---------- Feasible & Constraint Admissible Invariant Sets ---------- %%

d = size(sysd.B,2);n = size(sysd.B,1);
model = LTISystem('A', sysd.A+sysd.B*K,'Ts',T_sample);
% state constraints (A_x * x <= b_x)
X_set = Polyhedron([eye(n);-eye(n)],[xmax; -xmin]);
% input constraints (A_u * u <= b_u)
U_set = Polyhedron([eye(d); -eye(d)],[umax; -umin]);
% constraints admissible set
CA_set = Polyhedron([eye(d);-eye(d)]*K,[umax; -umin]);
% input constraint admissible set
IA_set = X_set&CA_set; % input constraint admissible set
% or use
% IA_set = Polyhedron(U_set.A*K_LMI,U_set.b)&X_set;
% invariant set with CA set
INVset_mpc = model.invariantSet('X',IA_set);
INVset_mpc_proj = projection(INVset_mpc,1:2);

MN_mpc = INVset_mpc.A; bN_mpc = INVset_mpc.b;
[W, Lu, c] = getWLc_TS(sysd.A, sysd.B, xmax, xmin, umax, umin, Gamma, Phi, MN_mpc, bN_mpc);

Px = Polyhedron('A',[-W Lu],'b',c);
Fset_x = projection(Px,1:2);             % Projection of Feasible set of X
%% -------------------- Closed Loop MPC Controller --------------------- %%

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
clear u_qp u_qp_list exitflag;

x = zeros(3,k_sim+1);
realX = x;
x(:,1) = x0;
realX(:,1) = x0+x_ss;

options_qp =  optimoptions('quadprog','Display','off');

for i = 1:k_sim
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
plot(0:k_sim,x_cloop','Linewidth',2);
grid on
set(gca,'FontWeight','bold')
% xlim([0 k_sim]);ylim([-2 3]);set(gca,'FontWeight','bold')
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',16,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold');ylabel('$States$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('State trajectory','FontSize',16,'FontWeight','bold');
subplot(2,1,2);
plot(0:k_sim-1,u_mpc_cloop,'Linewidth',2);
grid on
set(gca,'FontWeight','bold')
% % % xlim([0 k_sim]);ylim([-2 0]);set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16);ylabel('$u$','Interpreter','latex','FontSize',16);
title('Control Input','FontSize',16,'FontWeight','bold');

%% ------------------- Simulink MPC Controller Plot -------------------- %%
figure('Name','Simulink Results','NumberTitle','off');
subplot(2,1,1);
plot(0:k_sim,realX','Linewidth',2);
grid on
set(gca,'FontWeight','bold')
xlim([0 1000]);
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',16,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold'); ylabel('$States$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('State trajectory','FontSize',16,'FontWeight','bold');
subplot(2,1,2);
plot(0:k_sim-1,realU,'Linewidth',2);
grid on
xlim([0 1000]);
set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold');ylabel('$u$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('Control Input','FontSize',16,'FontWeight','bold');

%% -- Feasible and Invariant Constraint Admissible Set of States Plot -- %%
figure('Name','Feasible Set of x','NumberTitle','off');
set(gca,'FontWeight','bold')
plot(Fset_x); title('Trajectory Plot with Terminal Set','FontSize',16);
hold on
plot(INVset_mpc_proj,'color','yellow');
plot(realX(1,:),realX(2,:),'ogreen');
plot(x_cloop(1,:),x_cloop(2,:),'oblue');
xlabel('x_1','FontSize',14); ylabel('x_2','FontSize',14);
legend('Feasible Set','C.A. Invariant Set','Trajectory (Non-linear)','Trajectory (Linearized)','FontSize',14,'FontWeight','bold');
% figure("Name","Constraint Admissible Invariant Set",'NumberTitle','off')
% set(gca,'FontWeight','bold')
% plot(INVset_mpc_proj,'color','yellow'); title('Constraint Admissible Invariant Set of States','FontSize',14');
% xlabel('$x_1$','Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% 
% figure("Name","Feasible Set of States",'NumberTitle','off')
% set(gca,'FontWeight','bold')
% plot(Fset_x); title('Feasible Set of States','FontSize',14');
% xlabel('$x_1$','Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');