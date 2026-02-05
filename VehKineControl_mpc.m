clc;
close all;
clear all;
tic
% 车辆配置参数
VehConf.m  = 1573;
VehConf.Cf = 80000;
VehConf.Cr = 80000;
VehConf.Lf = 1.1;
VehConf.Lr = 1.58;
VehConf.Iz = 2873;

% 初始状态
x0    = 2;
y0    = 3;
phi0  = 0;

v0    = 2;
dT    = 0.01;
vofs1  = 0;
vofs2  = 0;

% 生成期望轨迹
[x1, y1, v1, phi1] = GenRefLineSegment(x0,      y0,      v0,      v0,     'left',    dT);  % 第一段路程
[x2, y2, v2, phi2] = GenRefLineSegment(x1(end), y1(end), v0,      v0+vofs1, 'left',  dT);  % 第二段路程
[x3, y3, v3, phi3] = GenRefLineSegment(x2(end), y2(end), v0+vofs2, v0,      'center', dT);  % 第三段路程

x_ref = [x1(1:end-1); x2(1:end-1); x3(1:end-1)];
y_ref = [y1(1:end-1); y2(1:end-1); y3(1:end-1)];
v_ref = [v1(1:end-1); v2(1:end-1); v3(1:end-1)];

phi_ref = [phi1(1:end-1); phi2(1:end-1); phi3(1:end-1)];

s_ref = cumsum([0; sqrt(diff(x_ref).^2 + diff(y_ref).^2)]);
len   = length(x_ref);
simT  = (0:(len-1))*dT;

% 初始状态空间变量
XCur = [y0; 0; phi0; 0];
VCur = v0;
Pos  = [x0; y0];

% 输入
delta   = 0;
acc     = 0;
Inputs  = [];
Outputs = [];
Ts      = [];

for i = 1:(length(simT)-10)
    
    t = simT(i);
    Inputs = [Inputs; delta, acc];
    Outputs = [Outputs; Pos(1), Pos(2), XCur(2), XCur(3), VCur];
    Ts = [Ts; t];
    
    % 根据距离查找最近点
    idx = find(s_ref > v0*t, 1);
    
    % 计算偏差
    x_err = Pos(1) - x_ref(idx);
    y_err = Pos(2) - y_ref(idx);
    phi = XCur(3);
    phi_err = phi - phi_ref(idx);
    
        % 道路曲率
    if idx < 5
        a = [x_ref(idx),  y_ref(idx)];
        b = [x_ref(idx+1),y_ref(idx+1)];
        c = [x_ref(idx+2),y_ref(idx+2)];
    else
        a = [x_ref(idx-1),y_ref(idx-1)];
        b = [x_ref(idx),  y_ref(idx)];
        c = [x_ref(idx+1),y_ref(idx+1)];
    end
    
    kDes = -getCurvature(a,b,c);  % 得到3个点曲率的函数
    
    ErrState = [x_err; y_err; phi_err];

    % 得到运动学MPC计算相关的A，B，C矩阵
    [Ad, Bd, Cd] = GetKineMPCControlMatrix(VehConf, VCur, phi_ref(idx), kDes*(VehConf.Lf+VehConf.Lr), dT);
    
       % 求解MPC问题
   % u = SolveLinearMPC(...
   %      A_k, B_k, C_k, ...
   %      diag([10,10,1]), diag([1,1]), ...
   %      [-1;-1], [1;1], ...
   %      ErrState, zeros(3*5,1), ...
   %      5);
   Dd = zeros(1,1);
%Ad3*3, Bd3*2, Cd3*2

% [Ad, Bd, Cd, Dd]=c2dm(Ac,Bc,Cc,Dc,dT,'zoh');
% tic
% Make the augmented state-space model
n1 = 0;
m1 = 0;
n_in = 0;

[m1, n1] = size(Cd);%3*2
[n1, n_in] = size(Bd);%3*2
A = eye(n1+m1, n1+m1);
A(1:n1, 1:n1) = Ad;
A(n1+1:n1+m1,1:n1) = Cd*Ad;
B = zeros(n1+m1, n_in);
B(1:n1,:) = Bd;
B(n1+1:n1+m1,:) = Cd*Bd;
C = zeros(m1, n1+m1);
C(:, n1+1:n1+m1) = eye(m1, m1);
D = Dd;


% Q=C'*C;
Q=1*eye(n1+m1,n1+m1);%5*5
R=1.0695*eye(n_in,n_in);%2*2

a1=0.55;
a2=0.55;
N1=3; %you can increase N
N2=4;
a11=[a1 a2];
N=[N1 N2];
Np=20;
Nc=10;

 % q=diag([10,10,1]);
 % r=diag([51,1]);
 % Q = [];
 % R = [];
 % 
 %    for i = 1:Np
 %        Q = blkdiag(Q, q); % 从三个不同大小的矩阵创建一个分块对角矩阵
 %        R = blkdiag(R, r); % 从三个不同大小的矩阵创建一个分块对角矩阵
 %    end


[M1,Lzerot]=Mdu(a11,N,n_in,Nc);
M0=Mu(a11,N,1,Nc);

u_min=-pi/6;
u_max=pi/6;
deltau_min=-pi/12;
deltau_max=-pi/12;
y_min=-100;
y_max=100;

u=0;

y_bar_min=y_min-y_ref;
y_bar_max=y_max-y_ref;





[Phi,F]=mpcgain(Ad,Bd,Cd,Nc,Np);
[Omega,Psi]=dmpc(A,B,a11,N,Np,Q,R);



%% constrains on U first
 % M=[M0;-M0;M1;-M1;Phi(1:Nc,1:2)*Lzerot;-Phi(1:Nc,1:2)*Lzerot];
 M=[M0;-M0;M1;-M1];

%closed-loop simulation performed recursively
%optimal control is also found recursively

% %% constrains on u first
% gamma=[(u_max-up)*ones(Nc,1);(-u_min+up)*ones(Nc,1);...
%     deltau_max*ones(Nc,1);-deltau_min*ones(Nc,1);...
%     y_bar_max(1:Nc,1)-F(1:Nc,1:4)*XCur;-y_bar_min(1:Nc,1)+F(1:Nc,1:4)*XCur];
gamma=[(u_max-u)*ones(Nc,1);(-u_min+u)*ones(Nc,1);...
    deltau_max*ones(Nc,1);-deltau_min*ones(Nc,1)];


Xf=[XCur; y_err];
% eta=QPhild(Omega,Psi(:,1:5)*Xf,M(1:40,:),gamma(1:40,:));
eta=QPhild(Omega,Psi(:,1:5)*Xf,M(1:40,:),gamma);

deltau=Lzerot*eta;

if (deltau>deltau_max)
    deltau=deltau_max;
end
if (deltau<deltau_min) 
    deltau=deltau_min;
end
u=u+deltau;
if (u>u_max)
    deltau=u_max-u; 
    u=u_max; 
end
if (u<u_min) 
    deltau=u_min-u; 
    u=u_min; 
end
    
    deltau = u(1) + kDes*(VehConf.Lf + VehConf.Lr);
    acc   = 0;
  
    [Pos, XCur, VCur] = GetNextPosition(VehConf, Pos, XCur, VCur, deltau, acc, dT);
    
end
toc
hold on;
plot(Outputs(:,1), Outputs(:,2), 'DisplayName', 'Actual Path')
plot(x_ref, y_ref,'DisplayName','Planning Path')
legend(gca,'show');
