function [ P ] = RotateAroundAxis( Q, axis, angle )
%RotateAroundAxis: rotate the geometry around defined axis with angle
cos_alpha=axis(3)/sqrt(axis(1)^2+axis(2)^2+axis(3)^2);
sin_alpha=sqrt(axis(1)^2+axis(2)^2)/sqrt(axis(1)^2+axis(2)^2+axis(3)^2);
cos_beta=axis(1)/sqrt(axis(1)^2+axis(2)^2);
sin_beta=axis(2)/sqrt(axis(1)^2+axis(2)^2);
if(axis(1)^2+axis(2)^2==0)
cos_beta=1;
sin_beta=0;
end

T1=[cos_alpha 0 -sin_alpha;
    0         1          0;
    sin_alpha 0  cos_alpha;];
T2=[cos_beta  sin_beta 0;
    -sin_beta cos_beta 0;
    0         0        1;];
T3=[cos(angle)  sin(angle) 0;
    -sin(angle) cos(angle) 0;
    0         0        1;];

P=Q*T2'*T1'*T3*T1*T2;

end

