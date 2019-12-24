function [p]=solver(A,B,C,F,e1,e2,e3,v)
%ABC是椭圆二次项参数，F是常数项参数
%e1,2,3是坐标初值
%v是车轮之间的空间坐标之差，如果没有默认为0，不进行筛选
%p是长度为5的向量，包含s1,s2,s3,R,t
    function [p]=f(x)
        alpha,d=14,0.01;
        %alpha为平面坐标与空间坐标变换参数
        % d为焦点到成像平面的距离
        % 二者均需要标定以后确定，修改函数中的参数
        s1=x(1);
        s2=x(2);
        s3=x(3);
        R=x(4);
        t=x(5);
        p=[s1 ^ 2 + s2 ^ 2 + s3 ^ 2 - 1, ...
            (e2 ^ 2 + e3 ^ 2 - R ^ 2) * s1 ^ 2 + (e2 * s2 + e3 * s3) ^ 2 - A * t / alpha ^ 2,...
            (e1 ^ 2 + e3 ^ 2 - R ^ 2) * s2 ^ 2 + (e1 * s1 + e3 * s3) ^ 2 - B * t / alpha ^ 2,...
            (e1 ^ 2 + e2 ^ 2 + e3 ^ 2 - R ^ 2) * s1 * s2 - (e1 * s1 + e2 * s2 + e3 * s3) * ...
            ( e1 * s2 + e2 * s1) - C * t / (2 * alpha ^ 2),...
            (e1 ^ 2 + e2 ^ 2 - R ^ 2) * s3 ^ 2 + (e1 * s1 + e2 * s2) ^ 2 - F * t / d ^ 2];
        
    end
%计算初值
factor = 1 / (A + C + (A ^ 2 + B ^ 2 + C ^ 2 + A * C) ^ 0.5);
begin_1 = (1 - C * factor) ^ 0.5;
begin_3 = B * factor / begin_1;
begin_2 = (1 - begin_1 ^ 2 - begin_3 ^ 2) ^ 0.5;
if v~=0
    p1 = abs(begin_1 * v(1) + begin_2 * v(2) + begin_3 * v(3));
    p2 = abs(begin_1 * v(1) - begin_2 * v(2) + begin_3 * v(3));
    if p1 > p2
        begin_2 = -begin_2;
    end
end
p=fsolve(f,[begin_1, begin_2, begin_3, 1, 2]);

end
        
       