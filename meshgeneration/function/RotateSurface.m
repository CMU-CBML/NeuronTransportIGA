function [ B ] = RotateSurface( A, norm_old, norm_new )
% Rotate the surface based on the norm vector change 
norm_old=norm_old/norm(norm_old);
norm_new=norm_new/norm(norm_new);
rotation_axis=cross(norm_old,norm_new);
alpha=acos(dot(norm_old,norm_new));
B=RotateAroundAxis(A,rotation_axis,alpha);
end

