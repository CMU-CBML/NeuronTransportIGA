function [ XYZ_smooth, D_smooth, tangent_vec] = BsplineSmooth( XYZ, D,seg_length, mode)
%Use the skeleton nodes as control points to construct B-spline
% and connect the points on the spline to create new series of segments.
% Input
% -----
% - XYZ: original skeleton nodes coordinates
% - D:   original skeleton nodes diameters
% - seg_length: approximate segment length for this branch, used to
%               calculate the sample point
% - mode: -- 1 for branch with both bifurcation ends
%         -- 2 for branch with one end as termination
%         -- 3 for the branch with start node
%
% Outputs
% -------
% - XYZ_smooth: smoothed skeleton nodes coordinates
% - D_smooth: smoothed skeleton nodes diameters
% - tangent_vec: the tangent vector at each node
%

[n_cpt,tmp]=size(XYZ);
d=0;

for i=2:n_cpt
    d=d+norm(XYZ(i,:)-XYZ(i-1,:));
end
u=zeros(n_cpt,1);
for k=2:n_cpt-1
    u(k)=u(k-1)+norm(XYZ(k,:)-XYZ(k-1,:))/d;
end

u(n_cpt)=1;


if (mode==1)
    sample_vec=[0 D(1)/d*1.5:seg_length/d:(1-1.0*D(end)/d) 1];
elseif mode==2
    sample_vec=[0 D(1)*1.5/d:seg_length/d:1];
    if(1-sample_vec(end)>0.05)
        sample_vec=[sample_vec 1];
    end
elseif mode==3
    sample_vec=[0:1.0*seg_length/d:(1-D(end)/d) 1];
elseif mode==4
    sample_vec=[0:1.0*seg_length/d:1]
    if(1-sample_vec(end)>0.05)
        sample_vec=[sample_vec 1];
    end
end
[tmp,n_sample]=size(sample_vec);

if n_cpt>3
    uu=u';
    knot=[0 0 0 0 uu(3:end-2) 1 1 1 1];
    if n_sample<4
        sample_vec=[0 0.45 0.67 1];
    end
elseif n_cpt==3
    knot=[0,0,0,1,1,1];
    if n_sample<4
        sample_vec=[0 0.45 0.67 1];
    end
else
    knot=[0,0,1,1];
    if n_sample<4
        sample_vec=[0 0.45 0.67 1];
    end
end

sp_xyz=spmak(knot,XYZ');
XYZ_smooth=fnval(sp_xyz,sample_vec);
XYZ_smooth=XYZ_smooth';

sp_d=spmak(knot,D');
D_smooth=fnval(sp_d,sample_vec);
D_smooth=D_smooth';

dsp_xyz=fnder(sp_xyz,1);
tangent_vec=fnval(dsp_xyz,sample_vec);
tangent_vec=tangent_vec';

end

