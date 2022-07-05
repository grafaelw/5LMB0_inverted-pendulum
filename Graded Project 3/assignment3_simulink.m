clear, clc,
%% -------------------- Simulation & Model Parameters ------------------ %%
%           Time
T_step =  1000;     % [s] 
pu = 0.5;           % [fr] step amplitude
T_sample = 0.004;   % [s]
T_sim = 1;          % []
k_sim = 5000;       % [] Quadprog looping simulation

R = 1;         %   [Ohm] Electrical resistance
L = 1/1000;    %   [H] Electrical inductance
k = 6/100;     %   [NA^-1] Motor constant
b = 1/1000;    %   [Nsm^-1] Friction coefficent
m = 7/100;     %   [kg] Pendulum mass 
l = 1/10;      %   [m] Pendulum length
J = m*l^2;     %   [kgm^2] Pendulum inertia
g = 9.81;      %   [ms^-2] Standard gravity

%% ----------------------- Initial Conditions -------------------------- %%
x1_0 = 1.1*pi/2; 
x2_0 = 10;   
x3_0 = -1.549;  

x0 = [x1_0;x2_0;x3_0];
x_ss = [1.133224715311108;0;-1.036669018180842];
u_ss = x_ss(3);
%% ---------------- Inactive Constraint MPC Controller ----------------- %%

% Inactive Constraints
N = 10; % Prediction Horizon
Q = eye(3); Q(1,1) = 120; Q(2,2) = 10; Q(3,3) = 10; % State weight
P = Q;
Re = 10;    % Control input weight

%% ----------------- Initial Condition MPC Controller ------------------ %%
B = [0;0;T_sample/L];

% Constrained MPC parameters
xmax = [2*pi;12;8];xmin = [-2*pi;-12;-8];
umax = 10;umin = -10;

%% --------------------- Simulink MPC Controller ---------------------- %%

x = zeros(3,k_sim+1);
z = zeros(3,k_sim+1);
% run("initialDataGenerator.m")
load("Uk-en-Rho.mat")

x(:,1) = x0;
z(:,1) = x0;

options_qp =  optimoptions('quadprog','Display','off');

for i = 1:k_sim
    warning off
    % P(k) matrix
    tic;
    % Phi and Gamma
    [Phi, Gamma] = LPVPhiGamma(B,N, Rho);

    % Prediction matrices for MPC
    [Psi, Omega] = QRPN2PsiOmega(Q,Re,P,N);
    G = 2*(Psi+Gamma'*Omega*Gamma);
    F = 2*Gamma'*Omega*Phi;
    [W, Lu, c] = getWLc(B,xmax,xmin,umax,umin,Gamma,Phi);
    %test: do this until U(i+1)-U(i)<epsillon
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
    u_input = u_mpc(i);% feeding the u_mpc to the simulink model
    sim('Pendulum_Nonlinear_System');
    x(:,i+1) = [x1_state(2,:); x2_state(2,:); x3_state(2,:)];
    z(:,i+1) = x(:,i+1)-x_ss;
    u_calibrated(i) = u_mpc(i)-u_ss;
    x1_0 = x1_state(2,:); x2_0 = x2_state(2,:); x3_0 = x3_state(2,:); 
    
    % calculating u and x
    
    predX(:,1) = x(:,i+1);

    for j = 2:N
        predX(:,j) = [1                     T_sample          0;
                      Rho(j)*m*g*l*T_sample/J (1-b*T_sample/J)  k*T_sample/J;
                       0,-k*T_sample/L,(1-R*T_sample/L)]*x(:,i)+B*u_qp(j);
    end
    
    % Calculating for the future P = {rho_i|k, ....., rho_N|k}
    Rho = sinc(x(:,i+1)/pi); 
    for ii=2:1:N-1
        Rho =  [Rho; sinc(predX(1,ii)/pi)];
    end
    timeLPVSimulink(i) = toc;
end
save("SimulinkLPV-time.mat",'timeLPVSimulink');
%% ------------------- Simulink MPC Controller Plot -------------------- %%
figure('Name','Simulink Results','NumberTitle','off');
subplot(2,1,1);
plot(0:k_sim,x','Linewidth',2);
grid on
set(gca,'FontWeight','bold')
%xlim([0 1000]);
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',16,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold'); ylabel('$States$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('State trajectory','FontSize',16,'FontWeight','bold');
subplot(2,1,2);
plot(0:k_sim-1,u_mpc,'Linewidth',2);
grid on
%xlim([0 1000]);
set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold');ylabel('$u$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('Control Input','FontSize',16,'FontWeight','bold');

%% ------------- Simulink MPC Controller Plot (Corrected) -------------- %%
figure('Name','Simulink Results Corrected','NumberTitle','off');
subplot(2,1,1);
plot(0:k_sim,z(1,:),'Linewidth',2);
grid on
hold on
plot(0:k_sim,z(2,:),'Linewidth',2);
plot(0:k_sim,z(3,:),'Linewidth',2);
set(gca,'FontWeight','bold')
% xlim([0 k_sim]);ylim([-2 3]);set(gca,'FontWeight','bold')
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',16,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold');ylabel('$States$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('State trajectory','FontSize',16,'FontWeight','bold');
subplot(2,1,2);
plot(0:k_sim-1,u_calibrated,'Linewidth',2);
grid on
set(gca,'FontWeight','bold')
% % % xlim([0 k_sim]);ylim([-2 0]);set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16);ylabel('$u$','Interpreter','latex','FontSize',16);
title('Control Input','FontSize',16,'FontWeight','bold');
