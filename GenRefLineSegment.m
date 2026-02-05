function [x, y, v, phi] = GenRefLineSegment(x0, y0, v0, v1, type, dT)
% 生成一段参考轨迹
% 参数:
% x0  起点的x坐标
% y0  起点的y坐标
% v0  起点的速度

% v1  终点的速度
% type 参考轨迹类型 [left | right | center]
% 如果类型是left，则生成的轨迹段是以当前点开始，向左边车道变道的轨迹
% 如果类型是right，则生成的轨迹段是以当前点开始，向右边车道变道的轨迹
% 如果类型是center，则生成的轨迹段是以当前点开始，直行的轨迹

% dT  生成轨迹点的时间间隔


% 输出
% x   生成轨迹点的x坐标
% y   生成轨迹点的y坐标
% v   生成轨迹点的速度
% phi 生成轨迹点的横摆角
    switch type
        case 'left'
            [x,y, v, phi] = GenRoadPoint(x0, y0, v0, v1, 4.0, dT, 3);% 真正生成轨迹的函数
        case 'right'
            [x,y, v, phi] = GenRoadPoint(x0, y0, v0, v1, -4.0, dT, 3);
        otherwise
            [x,y, v, phi] = GenRoadPoint(x0, y0, v0, v1, 0, dT, 3);
    end
end
function [x, y, v, phi] = GenRoadPoint(x0, y0, v0, v1, y1, dT, DurT)% 真正生成轨迹的函数
% 真正生成轨迹的函数
% 参数
% x0  初始点x坐标
% y0  初始点y坐标
% v0  初始点速度值
% v1  终点速度
% y1  终点y坐标
% dT  生成轨迹的时间间隔
% DurT 生成轨迹段的时间跨度
% 算法
% 采用横向5次纵向4次多项式拟合的算法
% 输出
% x   生成轨迹点的x坐标
% y   生成轨迹点的y坐标
% v   生成轨迹点的速度
% phi 生成轨迹点的横摆角
    cur_lateral_pos = 0;
    cur_lateral_spd = 0;
    cur_lateral_acc = 0;
    des_lateral_pos = y1;
    des_t = DurT;
    lat_qp = quintic_polynomial(cur_lateral_pos, cur_lateral_spd, cur_lateral_acc, des_lateral_pos, 0.0, 0.0, des_t);
    cur_lon_pos = 0;
    cur_speed   = v0;
    des_spd     = v1;
    lon_qp = quartic_polynomial(cur_lon_pos, cur_speed, 0.0, des_spd, 0.0, des_t);
    t = 0:dT:DurT;
    t = t';
    lat_qp = lat_qp(end:-1:1);
    lon_qp = lon_qp(end:-1:1);
    lon_pos = polyval(lon_qp, t);
    lon_spd = polyval([4,3,2,1].*lon_qp(1:end-1), t);
    lat_pos = polyval(lat_qp, t);
    lat_spd = polyval([5,4,3,2,1].*lat_qp(1:end-1),t);
    x = x0 + lon_pos;
    y = y0 + lat_pos;
    v = sqrt(lon_spd.^2 + lat_spd.^2);
    phi = atan(lat_spd./lon_spd);
end

function lat_qp = quintic_polynomial(xs, vxs, axs, xe, vxe, axe, T)
% 5次多项式拟合
    a0 = xs;
    a1 = vxs;
    a2 = axs/2;
    A = [T^3,   T^4,   T^5
         3*T^2, 4*T^3, 5*T^4
         6*T    12*T^2,20*T^3];
    b = [xe-a0-a1*T-a2*T^2
         vxe-a1-2*a2*T
         axe-2*a2];
    x = A\b;
    lat_qp = [a0, a1, a2, x(1),x(2),x(3)];
end


function lat_qp = quartic_polynomial(xs, vxs, axs, vxe, axe, T)
% 4次多项式拟合
    a0 = xs;
    a1 = vxs;
    a2 = axs/2;
    A = [3*T^2, 4*T^3;
         6*T    12*T^2,];
    b = [vxe-a1-2*a2*T
         axe-2*a2];
    x = A\b;
    lat_qp = [a0, a1, a2, x(1),x(2)];
end
