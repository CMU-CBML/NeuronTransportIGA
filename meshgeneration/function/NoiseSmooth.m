function [ B ] = NoiseSmooth( A, k, n )
%UseLaplace smoothing to remove noise nodes on the skeleton, both ends of
%the branch are fixed during the smoothing
% Input
% -----
% - A: old point location
% - k: smooth parameter
% - n: total iteration steps
% Outputs
% -------
% - B: new point location

tmp1=zeros(4,1);
tmp2=zeros(4,1);
delta=zeros(4,1);
[Ai Aj]=size(A);
for iterator=1:n
    for i=2:Ai-1
        tmp1=A(i,:);
        tmp2=0.5*(A(i-1,:)+A(i+1,:));
        delta=tmp2-tmp1;
        A(i,:)=tmp1+k*delta;
    end
end
B=A;
end

