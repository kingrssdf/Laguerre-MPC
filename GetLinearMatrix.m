function [A, B] = GetLinearMatrix(VehConf, V)
% 侧向状态空间模型中的A，B矩阵和车辆速度有关，这个函数的作用是根据车辆速度获取A，B矩阵

% 输入
% VehConf 车辆参数配置
% V 车辆速度


% 输出
% A A矩阵
% B B矩阵
    Cr = VehConf.Cr;
    Cf = VehConf.Cf;
    m  = VehConf.m;
    Lr = VehConf.Lr;
    Lf = VehConf.Lf;
    Iz = VehConf.Iz;
    
    
    A = zeros(4,4);
    A(1,2) = 1;
    
    A(2,2) = -2*(Cf+Cr)/(m*V);
    A(2,4) = -V - (2*Cf*Lf-2*Cr*Lr)/(m*V);
    
    A(3,4) = 1;
    
    A(4,2) = -(2*Cf*Lf-2*Cr*Lr)/(Iz*V);
    A(4,4) = -(2*Cf*Lf^2 + 2*Cr*Lr^2)/(Iz*V);
    
    B = zeros(4,1);
    B(2) = 2*Cf/m; 
    B(4) = 2*Lf*Cf/Iz;
end
