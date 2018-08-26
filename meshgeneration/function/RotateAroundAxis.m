function [ P ] = RotateAroundAxis( Q, axis, angle )
%RotateAroundAxis: rotate the geometry around defined axis with angle
axis=axis/norm(axis);
K=[0 -axis(3) axis(2);
   axis(3) 0 -axis(1);
   -axis(2) axis(1) 0];
I=eye(3);
R=I+sin(angle)*K+(1-cos(angle))*K*K;
P=Q*R';
end

