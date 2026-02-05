 


function [Phi, F] = mpcgain(Ap,Bp,Cp,Nc,Np);  %Tutorial 1.2
%function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e, Phi, F] = mpcgain(Ap,Bp,Cp,Nc,Np);
m1 = 0;
n1 = 0;
n_in = 0;
[m1, n1] = size(Cp);%3*3
[n1, n_in] = size(Bp);%3*2
% m11 = m1
% n11 = n1
% n_inn = n_in
A_e = eye(n1+m1, n1+m1);%6*6
A_e(1:n1, 1:n1) = Ap;%3*3
A_e(n1+1:n1+m1, 1:n1) = Cp*Ap;%3*3
B_e = zeros(n1+m1, n_in);%6*2
B_e(1:n1,:) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;
C_e = zeros(m1, n1+m1);%3*6
C_e(:,n1+1:n1+m1) = eye(m1,m1);
n = n1+m1;
h = zeros(Np,n1+m1);%9*6
F = zeros(Np,n1+m1);%5*6
h(1:3,:) = C_e;
F(1:3,:) = C_e*A_e;
for kk = 2:Np
    h(kk,:) = h(kk-1,:)*A_e;
    F(kk,:) = F(kk-1,:)*A_e;
end
v = h*B_e;
Phi = zeros(Np,Nc);  %declare the dimension of Phi
Phi(:,1:2) = v;        %first column of Phi
for i=2:Nc
    Phi(:,i) = [zeros(i-1,1); v(1:Np-i+1,1)];  %Toeplitz matrix
end
BarRs = ones(Np,1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;

%delta_u = inv(Phi_Phi + r_w*eye(Nc))*(Phi_R - Phi_F*x);
%Kmpc = [1 0 0..0_Nc]*inv(Phi_Phi + r_w*eye(Nc))*Phi_F;