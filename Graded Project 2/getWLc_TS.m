function [W, L, c] = getWLc_new( A, B, Xmax, Xmin, umax, umin, Gamma, Phi, M_N, b_N)
%GETWLC Summary of this function goes here
%   Detailed explanation goes here

nu = size(B,2);
nx = size(B,1);
N = size(Phi,1)/size(B,1);
           
Mi          = [zeros(nu,nx); 
               zeros(nu,nx);
               +eye(nx);
               -eye(nx)];
Ei          = [+eye(nu); -eye(nu); zeros(size(Mi,1)-2*nu,nu)];
bi          = [umax;
               -umin; 
               Xmax;
               -Xmin];

           
Dcal = [Mi;repmat(0*Mi,N-1,1);0*M_N];
Mcal = M_N;
for i = 2:N
    Mcal = blkdiag(Mi,Mcal);
end
Mcal = [zeros(size(Mi,1),size(Mcal,2));Mcal];
Ecal = Ei;
for i = 2:N
    Ecal = blkdiag(Ecal,Ei);
end
Ecal = [Ecal;zeros(size(M_N,1),size(Ecal,2))];
c = b_N;
for i = 1:N
    c = [bi;c];
end

L = Mcal*Gamma+Ecal;
W = -Dcal-Mcal*Phi;

end

