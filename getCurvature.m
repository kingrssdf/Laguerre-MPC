function curvature = getCurvature(a,b,c)
% 得到3个点曲率的函数
% 输入
% a a点坐标
% b b点坐标
% c c点坐标
    ab = b - a;
    bc = c - b;
    ac = c - a;
    sinC = (ac(2)*bc(1) - ac(1)*bc(2))/norm(ac)/norm(bc);
    curvature = 2*sinC/norm(ab);
end