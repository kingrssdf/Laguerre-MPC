function [A_k, B_k, C_k] = GetKineMPCControlMatrix(VehConf, vr, thetar, deltaru, ts)
% VehConf : 车辆参数
% vr      : 前向速度
% thetar  : 参考角度
% deltar  : 偏转角度
% ts      : 离散时间

     Lr = VehConf.Lr;
     Lf = VehConf.Lr;
    l = Lr + Lf;
    
    A_k = zeros(3, 3);
    A_k(1, 1) = 1;
    A_k(1, 3) = -vr*sin(thetar)*ts;
    A_k(2, 2) = 1;
    A_k(2, 3) = vr*cos(thetar)*ts;
    A_k(3, 3) = 1;
    
    B_k = zeros(3,2);
    B_k(1,1) = cos(thetar)*ts;
    B_k(2,1) = sin(thetar)*ts;
    B_k(3,1) = tan(deltaru)*ts/l;
    B_k(3,2) = vr*ts/(l*cos(deltaru)*cos(deltaru));
    
    C_k = eye(3,3);
end