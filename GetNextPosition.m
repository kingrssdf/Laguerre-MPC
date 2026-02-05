function [PosNext, XNext, VNext] = GetNextPosition(VehConf, Pos, XCur, VCur, delta, acc, dT)
% 根据当前车辆状态，计算下一个状态

% 输入
% VehConf 车辆配置参数
% Pos     车辆位置
% XCur    状态变量
% VCur    当前速度
% delta   转向角输入
% acc     纵向加速度
% dT      时间步长


% 输出
% PosNext  下一时刻的位置
% XNext    下一时刻的状态变量
% VNext    下一时刻的速度

% N是分段值
N = 1;
dt = dT/N;

v = VCur;
x = XCur;

xpos = Pos(1);
ypos = Pos(2);

for i = 1:N
    phi = x(3);
    xpos = xpos + v*dt*cos(phi);
    ypos = ypos + v*dt*sin(phi);
    
    vy = x(2);
%     xpos = xpos + vy*dt*cos(phi-pi/2);
%     ypos = ypos + vy*dt*sin(phi-pi/2);
    
    v = v + acc*dT/N;
    [A, B] = GetLinearMatrix(VehConf, v);
    
    % 这里面包含离散化的过程
    x = (A*dt+eye(4))*x + B*dt*delta;
end

XNext = x;
VNext = v;
PosNext = [xpos; ypos];
end
