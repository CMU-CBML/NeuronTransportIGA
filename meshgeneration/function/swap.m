function [ B ] = swap( A, i, j )

B=A;
B(i,:)=A(j,:);
B(j,:)=A(i,:);

end

