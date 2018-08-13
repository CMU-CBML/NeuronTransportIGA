function [ B ] = swap( A, i, j )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
B=A;
B(i,:)=A(j,:);
B(j,:)=A(i,:);

end

