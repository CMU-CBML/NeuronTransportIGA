function [ Q ] = Projection( A, B, C, P )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n_vec=cross(B-A,C-A)/norm(cross(B-A,C-A));
a=n_vec(1);
b=n_vec(2);
c=n_vec(3);
d=-(a*A(1)+b*A(2)+c*A(3));

k=(a*P(1)+b*P(2)+c*P(3)+d)/(a^2+b^2+c^2);
Q(1)=P(1)-k*a;
Q(2)=P(2)-k*b;
Q(3)=P(3)-k*c;
end

