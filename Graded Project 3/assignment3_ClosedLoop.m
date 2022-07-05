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

x0 = [0.55*pi;10;-1.549];
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
xmax = [2*pi;12;8];
xmin = [-2*pi;-12;-8];
umax = 10;umin = -10;

%% -------------------- Closed Loop MPC Controller --------------------- %%
% x_ss = [theta;0;-(m*g*l*sin(theta))/k;];    % Steady States of X
% u_ss = x_ss(3);             % Steady State of U
x = zeros(3,k_sim+1);
z = zeros(3,k_sim+1);
% run("initialDataGenerator.m")
load("Uk-en-Rho.mat")
Rho0 = Rho; % from called function of Rho&predX
% P, Q, R, N, x(k),{x(2|k-1),....,x(N|k-1)} are given
% x(0|k) = x0
x(:,1) = x0;
z(:,1) = x0;
% x(1|k) = A(rho(0|k)*x(0|k)+B*u(0|k)  B is rho independent
% rho(0|k) = fbar(x(k));
% rho_zero = sinc(x0(1)/pi);
% Rho = rho_zero;
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
    fprintf("Closed-loop Iteration: %d \n",i);
%     if norm(u_qp-Uk,1) <= 0.01
%         break
%     else
%         Uk = u_qp;
%     end
    u_qp_list(:,i) = u_qp;
    v_mpc(i) = u_qp(1);
    % calculating u and x

     x(:,i+1) = [1                     T_sample          0;
     Rho(1)*m*g*l*T_sample/J (1-b*T_sample/L)  k*T_sample/J;
     0                   -k*T_sample/L  (1-R*T_sample/L)]*x(:,i) + B*v_mpc(i); % Closed-loop MPC

     z(:,i+1) = x(:,i+1)-[1                     T_sample          0;
     Rho(1)*m*g*l*T_sample/J (1-b*T_sample/L)  k*T_sample/J;
     0                   -k*T_sample/L  (1-R*T_sample/L)]*x_ss - B*u_ss; % Closed-loop corrected MPC

    corrected_u(i) = v_mpc(i)-u_ss;
    
    predX(:,1) = x(:,i+1);

    for j = 2:N
        predX(:,j) = [1                     T_sample          0;
                      Rho(j)*m*g*l*T_sample/J (1-b*T_sample/J)  k*T_sample/J;
                       0,-k*T_sample/L,(1-R*T_sample/L)]*x(:,i)+B*u_qp(j);
    end
    
    Rho = sinc(x(:,i+1)/pi);

    for ii=2:1:N-1
        Rho =  [Rho; sinc(predX(1,ii)/pi)];
    end
    

    timeLPV_CL(i) = toc;

%     if norm(Rho0-Rho,1) <= 0.01
%         Rho = Rho;
%     else
%         
%     end
    
%     u_mpc_cloop(i) = v_mpc(i) + u_ss;
%     x_cloop(:,i+1) = z(:,i+1) + sysd.A*x_ss+ sysd.B*u_ss;
end

save("ClosedLoopLPV-time.mat",'timeLPV_CL');

%% ------------------ Closed-loop MPC Controller Plot ------------------ %%
figure('Name','Closed Loop MPC','NumberTitle','off');
subplot(2,1,1);
plot(0:k_sim,x(1,:),'Linewidth',2);
grid on
hold on
plot(0:k_sim,x(2,:),'Linewidth',2);
plot(0:k_sim,x(3,:),'Linewidth',2);
set(gca,'FontWeight','bold')
% xlim([0 k_sim]);ylim([-2 3]);set(gca,'FontWeight','bold')
legend('$x_1$','$x_2$','$x_3$','Interpreter','latex','FontSize',16,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16,'FontWeight','bold');ylabel('$States$','Interpreter','latex','FontSize',16,'FontWeight','bold');
title('State trajectory','FontSize',16,'FontWeight','bold');
subplot(2,1,2);
plot(0:k_sim-1,v_mpc,'Linewidth',2);
grid on
set(gca,'FontWeight','bold')
% % % xlim([0 k_sim]);ylim([-2 0]);set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16);ylabel('$u$','Interpreter','latex','FontSize',16);
title('Control Input','FontSize',16,'FontWeight','bold');


%% ------------------ Closed-loop MPC Controller Plot (Corrected) ------------------ %%
figure('Name','Closed Loop MPC Corrected','NumberTitle','off');
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
plot(0:k_sim-1,corrected_u,'Linewidth',2);
grid on
set(gca,'FontWeight','bold')
% % % xlim([0 k_sim]);ylim([-2 0]);set(gca,'FontWeight','bold')
xlabel('$k$','Interpreter','latex','FontSize',16);ylabel('$u$','Interpreter','latex','FontSize',16);
title('Control Input','FontSize',16,'FontWeight','bold');
