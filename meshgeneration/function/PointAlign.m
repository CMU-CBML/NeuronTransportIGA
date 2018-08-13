function [ D ] = PointAlign( A, B, C )

k=dot((A-B),(C-B))/norm(A-B)^2;
D=k*A+(1-k)*B;
D=(A+B)/2;
end

