clc;
clear;
addpath(genpath(pwd));
start_trees;
%% User settings


% Input and output path setting
% io_path='..//example//cylinder//'; % Figure 3E
io_path='..//example//bifurcation//'; % Figure 3F
% io_path='..//example//3bifurcation//'; % Figure 3G
% io_path='..//example//movie2//'; % Figure 5
% io_path='..//example//movie5//'; % Figure 6
% io_path='..//example//cell3traceRN1//';% Figure 7A
% io_path='..//example//nelson2//';% Figure 7D
% io_path='..//example//purkinje//';% Figure 7G

parameter_file=[io_path,'mesh_parameter.txt'];
skeleton_input=[io_path,'skeleton_smooth.swc'];
velocity_output=[io_path,'initial_velocityfield.txt'];
hex_output=[io_path,'controlmesh.vtk'];

% Read skeleton information
load_tree(skeleton_input);

% parameter setting
var=LoadParameter(parameter_file);
ratio_refine=var(5); % the parameter used to calculate the refinement around bifurcation region
%% Extract skeleton information and initialize labels
location=[trees{1}.X,trees{1}.Y,trees{1}.Z];
d=trees{1}.D;
id=idpar_tree';
id(1)=0;
ipar=ipar_tree();

branch=B_tree();
termination=T_tree();

[nx,ny]=size(id);
[m,n]=size(find(trees{1}.dA~=0));
[mm,nn]=size(ipar);

in_point_label=-1;
wall_label=0;
tip_label=1;
branch_label=-2;

ijk_label=zeros(nx,1); % i-0 j-1 k-2

%% Read template mesh
template_p=load('.//template//template_circle_points90.txt');
template_half_p=load('.//template//template_halfcircle_points90.txt');
template_merge_p=load('.//template//template_merge120_points90.txt');

template_e=load('.//template//template_circle_elements.txt');
template_half_e=load('.//template//template_halfcircle_elements.txt');
template_merge_e=load('.//template//template_merge_elements.txt');
template_merge_e_reverse=template_merge_e;

[m1,n1]=size(template_p); %template points number
[m2,n2]=size(template_e); %template elements number
[m3,n3]=size(template_merge_p); %branch template points number
[m4,n4]=size(template_merge_e); %branch template elements number

for i=1:m1
    velocity_value(i)=(1-norm(template_p(i,:))^2);
end
inner_point_index=find(sqrt(template_p(:,1).^2+template_p(:,2).^2+template_p(:,3).^2)<0.95);
boundary_point_index=find(sqrt(template_p(:,1).^2+template_p(:,2).^2+template_p(:,3).^2)>0.95);
boundary_merge_point_index=find(sqrt(template_merge_p(:,1).^2+template_merge_p(:,2).^2+template_merge_p(:,3).^2)>0.95);

for i=1:m4
    if i>m4/3*2
        template_merge_e_reverse(i,:)=[template_merge_e(i,1),template_merge_e(i,2),template_merge_e(i,5),template_merge_e(i,4),template_merge_e(i,3)];
    end
end

%% Construct bifurcation template
tmp_merge_e_order=template_merge_e_reverse;
% modify branch element order
a=[90 91 92 93 95 96 97 100 101 105 125:134];
b=[114 109 104 99 113 108 103 112 107 111 124:-2:116 123:-2:115];
a=a+1+m2/2;
b=b+1+m2/2;
aa=a+m2/4;
bb=b+m2/4;
[m,n]=size(a);
for i=1:n
    tmp_merge_e_order=swap(tmp_merge_e_order,a(i),b(i));
    tmp_merge_e_order=swap(tmp_merge_e_order,aa(i),bb(i));
end
for i=1:m2/4
    tmp_merge_e_order=swap(tmp_merge_e_order,i+m2,i+m2*5/4);
end

% determine the boundary point order in template
botbc_index=[];
leftbc_index=[];
rightbc_index=[];

for i=1:m3
    if(norm(template_merge_p(i,:))>0.95)
        if(template_merge_p(i,2)>0)
            leftbc_index=[leftbc_index;i];
        elseif(template_merge_p(i,2)<0)
            rightbc_index=[rightbc_index;i];
        end
        if(template_merge_p(i,3)>0)
            botbc_index=[botbc_index;i];
        end
    end
end

extra_pt=[2 25 124 217;
    4 32 125 218]; % extraordinary points

% assign template element order to branch
template_branch_bottom=template_e;
template_branch_left=[template_e(1:m2/2,:);template_merge_e((m2+1):end,:)];
% template_branch_right=[template_merge_e((m2+1):end,:);template_e(m4+1:end,:)];
template_branch_right=[tmp_merge_e_order((m2+1):end,:);template_e(m2/2+1:end,:)];

template_branch_pt_bot=template_merge_p(1:m1,:);
template_branch_pt_left=template_merge_p(1:m1,:);
template_branch_pt_right=template_merge_p(1:m1,:);
template_branch_pt_bot_terminal=template_p;

rotate90_cw120=[1 0 0;0 cosd(-120) -sind(-120);0 sind(-120) cosd(-120)];
rotate90_ccw120=[1 0 0;0 cosd(120) -sind(120);0 sind(120) cosd(120)];

for i=1:m1
    if template_merge_p(i,2)<0
        template_branch_pt_left(i,:)=template_merge_p(i,:)*rotate90_ccw120;
    end
    if template_merge_p(i,2)>0
        template_branch_pt_right(i,:)=template_merge_p(i,:)*rotate90_cw120;
    end
end

%% Define the temporary point and element in each layer
tmp_p=template_p;
tmp_mp=template_merge_p;
tmp_e=zeros(m2,8);
tmp_me=zeros(m4,8);
tmp_branch_be=zeros(m2,8);
tmp_branch_le=zeros(m2,8);
tmp_branch_re=zeros(m2,8);
n_layer=zeros(ny,ny);
branch_element_location=cell(ny);
Segment_Vector=cell(ny,ny);
Tangent_Vector=cell(ny,3);
NodeLayer=cell(ny,4);
BranchLayerPoint=cell(ny,3);

AllPoint=[];
LayerElement=[];
AllElement=[];
AllLabel=[];
AllVelocity=[];
AllLayerIndex=[];

for i=1:ny
    for j=1:ny
        if trees{1}.dA(i,j)~=0
            Segment_Vector{i,j}=location(i,:)-location(j,:);
            sv=location(i,:)-location(j,:);
            %             n_layer(i,j)=ceil(norm(Segment_Vector{i,j})/(d(j)*1.1));
            %             if(termination(i) || j==1)
            %                 n_layer(i,j)=ceil(norm(Segment_Vector{i,j})/(d(j)*0.25));
            %             end
            % set layer number for normal segments
            n_layer(i,j)=1;
            % set layer number segments around bifurcation for refinement
            % the bifurcation diameter is used to calculate the total layer numbers
            if(branch(i))
                n_layer(i,j)=ceil(norm(Segment_Vector{i,j})/(d(i)*ratio_refine));
            end
            if(branch(j))
                n_layer(i,j)=ceil(norm(Segment_Vector{i,j})/(d(j)*ratio_refine));
            end
        end
    end
end

%% Skeleton segmentation
% Group the all branches to two catergories:
% -- 1. branch between two bifurcation points;
% -- 2. branch with one bifurcation point and one termination end
[sect, vec]=dissect_tree(trees{1});
[n_sect, tmp]= size(sect);
[n_bif,tmp]=size(find(branch==1));
sect_pt=cell(n_sect,1); %save pt index for each section (may include branch node)
bif_term_pt=cell(n_sect,1); % save pt index for bif-term section (branch node excluded)
% bif_pt=cell(n_bif,1); %save pt index for bif region,
bif_pt=zeros(n_bif,4);
bif_ele=zeros(n_bif,4);%save ele index for bif region,
%order:[bifur_node, par_node, bifur1_node, bifur2_node]

for i=1:n_sect
    tmp_par=ipar(sect(i,2),:);
    loc=find(tmp_par==sect(i,1));
    sect_pt{i}=tmp_par(loc:-1:1);
end

index_bif=1;
for index_sec=1:n_sect
    [sect_ptnum,tmp]=size(sect_pt{index_sec}');
    index_start_pt=sect_pt{index_sec}(1);
    index_end_pt=sect_pt{index_sec}(sect_ptnum);
    
    bif_term_start=1;
    bif_term_end=sect_ptnum;
    
    if(branch(index_end_pt))
        par_node=id(index_end_pt);
        ans_tmp=find(trees{1}.dA(:,index_end_pt)==1);
        bifur1_node=ans_tmp(1);
        bifur2_node=ans_tmp(2);
        bif_pt(index_bif,:)=[index_end_pt,par_node,bifur1_node,bifur2_node];
        index_bif=index_bif+1;
        bif_term_end=sect_ptnum-1;
    end
    if(branch(index_start_pt))
        bif_term_start=2;
    end
    bif_term_pt{index_sec}=sect_pt{index_sec}( bif_term_start: bif_term_end);
end

%% Generate layers for all bifurcation regions
for index_bif=1:n_bif
    i=bif_pt(index_bif,1);
    par_node=bif_pt(index_bif,2);
    bifur1_node=bif_pt(index_bif,3);
    bifur2_node=bif_pt(index_bif,4);
    sv_parpar=Segment_Vector{par_node,id(par_node)};
    sv_par=Segment_Vector{i,par_node};ri=d(par_node)/2.0;
    sv_bifur1=Segment_Vector{bifur1_node,i};rj=d(bifur1_node)/2.0; %rj=(ri+rj)/2;
    sv_bifur2=Segment_Vector{bifur2_node,i};rk=d(bifur2_node)/2.0; %rk=(ri+rk)/2;
    ijk_label(bifur1_node)=1;
    ijk_label(bifur2_node)=2;
    
    tmp_mp=template_merge_p*(ri+rj+rk)/3;
    tmp_bp=template_branch_pt_bot* (ri+rj+rk)/3;
    tmp_lp=template_branch_pt_left*(ri+rj+rk)/3;
    tmp_rp=template_branch_pt_right*(ri+rj+rk)/3;
    
    
    tmp_bpt=template_branch_pt_bot_terminal*ri;
    tmp_lpt=template_branch_pt_bot_terminal*rj;
    tmp_rpt=template_branch_pt_bot_terminal*rk;
    
    %Using local coordinates to calculate three half planes
    %Bifurcation point is origin ([0,0,0])
    vi=-sv_par/norm(sv_par);
    vj=sv_bifur1/norm(sv_bifur1);
    vk=sv_bifur2/norm(sv_bifur2);
    
    aij=acos(dot(vi,vj)/norm(vi)/norm(vj));
    aik=acos(dot(vi,vk)/norm(vi)/norm(vk));
    akj=acos(dot(vk,vj)/norm(vk)/norm(vj));
    
    Kij=(ri*vj+rj*vi)/norm((ri*vj+rj*vi));
    Kik=(ri*vk+rk*vi)/norm((ri*vk+rk*vi));
    Kkj=(rk*vj+rj*vk)/norm((rk*vj+rj*vk));
    
    if(dot(cross(vi, Kij),cross(vi,Kik))>0)
        if aij>aik
            Kij=-Kij;
        else
            Kik=-Kik;
        end
    end
    
    if(aij<=pi/2)
        spij=Kij*ri/sin(atan(ri/rj));
    else
        spij=Kij*(ri+rj)/2;
    end
    if(aik<=pi/2)
        spik=Kik*ri/sin(atan(ri/rk));
    else
        spik=Kik*(ri+rk)/2;
    end
    if(akj<=pi/2)
        spkj=Kkj*rk/sin(atan(rk/rj));
    else
        spkj=Kkj*(rk+rj)/2;
    end
    cpn=cross((spik-spkj),(spij-spkj))/norm(cross((spik-spkj),(spij-spkj)));
    fprintf('i=%d original_cpn=[%f %f %f]\n',i,cpn(1),cpn(2),cpn(3));
    
    % Branch Template rotation
    cp1=cpn*(ri+rj+rk)/3;
    cp2=-cpn*(ri+rj+rk)/3;
    
    fprintf('i=%d cpn=[%f %f %f]\n',i,cpn(1),cpn(2),cpn(3));
    
    w_bpt=cross(sv_par,cp1);w_bpt=w_bpt/norm(w_bpt);
    w_lpt=cross(sv_bifur1,cp1);w_lpt=w_lpt/norm(w_lpt);
    w_rpt=cross(sv_bifur2,cp1);w_rpt=w_rpt/norm(w_rpt);
    for k=1:m1
        tmp_bpt(k,:)=cp1/norm(cp1)*tmp_bpt(k,1)+w_bpt*tmp_bpt(k,2);
        tmp_lpt(k,:)=cp1/norm(cp1)*tmp_lpt(k,1)+w_lpt*tmp_lpt(k,2);
        tmp_rpt(k,:)=cp1/norm(cp1)*tmp_rpt(k,1)+w_rpt*tmp_rpt(k,2);
    end
    
    
    %Deal with half plane separate points
    ciji=w_bpt*ri-sv_par;
    ciki=-w_bpt*ri-sv_par;
    cijj=w_lpt*rj+sv_bifur1;
    cjkj=-w_lpt*rj+sv_bifur1;
    cjkk=w_rpt*rk+sv_bifur2;
    cikk=-w_rpt*rk+sv_bifur2;
    
    ni_layer=n_layer(i,par_node);
    nj_layer=n_layer(bifur1_node,i);
    nk_layer=n_layer(bifur2_node,i);
    
    ciji=(spij*(ni_layer-1)+ciji)/ni_layer;
    ciki=(spik*(ni_layer-1)+ciki)/ni_layer;
    cijj=(spij*(nj_layer-1)+cijj)/nj_layer;
    cjkj=(spkj*(nj_layer-1)+cjkj)/nj_layer;
    cikk=(spik*(nk_layer-1)+cikk)/nk_layer;
    cjkk=(spkj*(nk_layer-1)+cjkk)/nk_layer;
    
    kij1=sqrt(norm(ciji)^2+norm(spij)^2-2*dot(spij,ciji));
    kij2=sqrt(norm(cijj)^2+norm(spij)^2-2*dot(spij,cijj));
    spij=(ciji*kij2+cijj*kij1)/(kij1+kij2);
    
    kik1=sqrt(norm(ciki)^2+norm(spik)^2-2*dot(spik,ciki));
    kik2=sqrt(norm(cikk)^2+norm(spik)^2-2*dot(spik,cikk));
    spik=(ciki*kik2+cikk*kij1)/(kik1+kik2);
    
    kjk1=sqrt(norm(cjkj)^2+norm(spkj)^2-2*dot(spkj,cjkj));
    kjk2=sqrt(norm(cjkk)^2+norm(spkj)^2-2*dot(spkj,cjkk));
    spkj=(cjkj*kjk2+cjkk*kjk1)/(kjk1+kjk2);
    
    for k=1:m1
        if(template_branch_pt_bot(k,2)<0)
            tmp_bp(k,:)=cp1*template_branch_pt_bot(k,1)+spik*norm(template_branch_pt_bot(k,2:3));
        elseif(template_branch_pt_bot(k,2)>0)
            tmp_bp(k,:)=cp1*template_branch_pt_bot(k,1)+spij*norm(template_branch_pt_bot(k,2:3));
        else
            tmp_bp(k,:)=cp1*template_branch_pt_bot(k,1);
        end
        if(template_branch_pt_left(k,3)<0)
            tmp_lp(k,:)=cp1*template_branch_pt_left(k,1)+spij*norm(template_branch_pt_left(k,2:3));
        elseif(template_branch_pt_left(k,3)>0)
            tmp_lp(k,:)=cp1*template_branch_pt_left(k,1)+spkj*template_branch_pt_left(k,3);
        else
            tmp_lp(k,:)=cp1*template_branch_pt_left(k,1);
        end
        if(template_branch_pt_right(k,3)<0)
            tmp_rp(k,:)=cp1*template_branch_pt_right(k,1)+spik*norm(template_branch_pt_right(k,2:3));
        elseif(template_branch_pt_right(k,3)>0)
            tmp_rp(k,:)=cp1*template_branch_pt_right(k,1)+spkj*template_branch_pt_right(k,3);
        else
            tmp_rp(k,:)=cp1*template_branch_pt_right(k,1);
        end
    end
    
    for k=1:m3
        if(template_merge_p(k,3)>0) %plane separate jk
            tmp_mp(k,:)=cp1*template_merge_p(k,1)+spkj*template_merge_p(k,3);
        elseif template_merge_p(k,2)>0 %plane separate ij
            tmp_mp(k,:)=cp1*template_merge_p(k,1)+spij*norm(template_merge_p(k,2:3));
        elseif template_merge_p(k,2)<0 %plane separate ik
            tmp_mp(k,:)=cp1*template_merge_p(k,1)+spik*norm(template_merge_p(k,2:3));
        else
            tmp_mp(k,:)=cp1*template_merge_p(k,1);
        end
    end
    
    
    % Move points to branch location
    for inode=1:m3
        tmp_mp(inode,:)= tmp_mp(inode,:)+location(i,:)-location(1,:);
    end
    for inode=1:m1
        tmp_bp(inode,:)= tmp_bp(inode,:)+location(i,:)-location(1,:);
        tmp_lp(inode,:)= tmp_lp(inode,:)+location(i,:)-location(1,:);
        tmp_rp(inode,:)= tmp_rp(inode,:)+location(i,:)-location(1,:);
        tmp_bpt(inode,:)=tmp_bpt(inode,:)+location(i,:)-location(1,:)-sv_par;
        tmp_lpt(inode,:)=tmp_lpt(inode,:)+location(i,:)-location(1,:)+sv_bifur1;
        tmp_rpt(inode,:)=tmp_rpt(inode,:)+location(i,:)-location(1,:)+sv_bifur2;
    end
    
    % Save Branch points and Branch layer element
    tmp_label=ones(m3,1);
    tmp_label=in_point_label*tmp_label;
    tmp_label(boundary_merge_point_index,:)=wall_label;
    for ii=1:m3
        tmp_veloctity_branch(ii,1:3)=[0 0 0];
    end
    
    [tmp_pointnumber,tmp]=size(AllPoint);
    bif_ele(index_bif,1)=tmp_pointnumber;% record the element information around bifurcation
    AllPoint=[AllPoint;tmp_mp];
    AllLabel=[AllLabel;tmp_label];
    AllVelocity=[AllVelocity;tmp_veloctity_branch];
    
    for k=1:m4
        tmp_me(k,1:4)=template_merge_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
    end
    for k=1:m2
        tmp_branch_be(k,1:4)=template_branch_bottom(k,2:5)+tmp_pointnumber*[1,1,1,1];
        tmp_branch_le(k,1:4)=template_branch_left(k,2:5)+tmp_pointnumber*[1,1,1,1];
        tmp_branch_re(k,1:4)=template_branch_right(k,2:5)+tmp_pointnumber*[1,1,1,1];
    end
    LayerElement=[LayerElement;tmp_me(:,1:4)];
    NodeLayer{i,1}=tmp_mp;
    NodeLayer{i,2}=tmp_branch_be(:,1:4);
    NodeLayer{i,3}=tmp_branch_le(:,1:4);
    NodeLayer{i,4}=tmp_branch_re(:,1:4);
    BranchLayerPoint{i,1}=tmp_bp;
    BranchLayerPoint{i,2}=tmp_lp;
    BranchLayerPoint{i,3}=tmp_rp;
    
    tmp_label=ones(m1,1);
    tmp_label=in_point_label*tmp_label;
    tmp_label(boundary_point_index,:)=wall_label;
    
    % Save bottom terminal points and element
    if(par_node==1)
        tmp_label=ones(m1,1);
        tmp_label=tip_label*tmp_label;
        tmp_label(boundary_point_index,:)=wall_label;
        tip_label=tip_label+1;
    end
    for ii=1:m1
        tmp_veloctity(ii,1:3)=sv_parpar/norm(sv_parpar)*norm(velocity_value(ii));
    end
    [tmp_pointnumber,tmp]=size(AllPoint);
    AllPoint=[AllPoint;tmp_bpt];
    AllLabel=[AllLabel;tmp_label];
    AllVelocity=[AllVelocity; tmp_veloctity];
    for k=1:m2
        tmp_e(k,1:4)=template_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
    end
    LayerElement=[LayerElement;tmp_e(:,1:4)];
    NodeLayer{par_node,1}=tmp_bpt;
    NodeLayer{par_node,2}=tmp_e(:,1:4);
    NodeLayer{par_node,3}=tmp_e(:,1:4);
    NodeLayer{par_node,4}=cp1;
    % record the element information around bifurcation
    bif_ele(index_bif,2)=tmp_pointnumber;
    
    
    % Save left terminal points and element
    if(termination(bifur1_node))
        tmp_label=ones(m1,1);
        tmp_label=tip_label*tmp_label;
        tmp_label(boundary_point_index,:)=wall_label;
        tip_label=tip_label+1;
        
    end
    for ii=1:m1
        tmp_veloctity(ii,1:3)=sv_bifur1/norm(sv_bifur1)*norm(velocity_value(ii));
    end
    [tmp_pointnumber,tmp]=size(AllPoint);
    AllPoint=[AllPoint;tmp_lpt];
    AllLabel=[AllLabel;tmp_label];
    AllVelocity=[AllVelocity; tmp_veloctity];
    for k=1:m2
        tmp_e(k,1:4)=template_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
    end
    LayerElement=[LayerElement;tmp_e(:,1:4)];
    NodeLayer{bifur1_node,1}=tmp_lpt;
    NodeLayer{bifur1_node,2}=tmp_e(:,1:4);
    NodeLayer{bifur1_node,3}=tmp_e(:,1:4);
    NodeLayer{bifur1_node,4}=cp1;
    % record the element information around bifurcation
    bif_ele(index_bif,3)=tmp_pointnumber;
    
    % Save right terminal points and element
    if(termination(bifur2_node))
        tmp_label=ones(m1,1);
        tmp_label=tip_label*tmp_label;
        tmp_label(boundary_point_index,:)=wall_label;
        tip_label=tip_label+1;
    end
    for ii=1:m1
        tmp_veloctity(ii,1:3)=sv_bifur2/norm(sv_bifur2)*norm(velocity_value(ii));
    end
    [tmp_pointnumber,tmp]=size(AllPoint);
    AllPoint=[AllPoint;tmp_rpt];
    AllLabel=[AllLabel;tmp_label];
    AllVelocity=[AllVelocity; tmp_veloctity];
    for k=1:m2
        tmp_e(k,1:4)=template_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
    end
    LayerElement=[LayerElement;tmp_e(:,1:4)];
    NodeLayer{bifur2_node,1}=tmp_rpt;
    NodeLayer{bifur2_node,2}=tmp_e(:,1:4);
    NodeLayer{bifur2_node,3}=tmp_e(:,1:4);
    NodeLayer{bifur2_node,4}=cp1;
    % record the element information around bifurcation
    bif_ele(index_bif,4)=tmp_pointnumber;
end

%% Calculate layers for bif-term branches and bif-bif branches
if n_bif~=0
    for index_sec=1:n_sect
        [sect_ptnum,tmp]=size(bif_term_pt{index_sec}');
        index_start_pt=bif_term_pt{index_sec}(1);
        index_end_pt=bif_term_pt{index_sec}(sect_ptnum);
        n_insert_layer=sect_ptnum-2;
        % bif-term segments
        if(index_start_pt==1)
            n_insert_layer=n_insert_layer+1;
            sec_vec_end=Segment_Vector{index_end_pt,id(index_end_pt)};
            ref_vec_end=NodeLayer{index_end_pt,4};
            w=cross(sec_vec_end,ref_vec_end); w=w/norm(w);
            ref_vec_next=ref_vec_end/norm(ref_vec_end);
            for index_sec_pt=sect_ptnum:-1:2
                
                j=bif_term_pt{index_sec}(index_sec_pt-1);
                i=bif_term_pt{index_sec}(index_sec_pt);
                sv=Segment_Vector{i,j};
                ref_vec_next=cross(w,sv); ref_vec_next=ref_vec_next/norm(ref_vec_next);
                
                for ii=1:m1
                    tmp_p(ii,:)=template_p(ii,1)*ref_vec_next+template_p(ii,2)*w;
                end
                tmp_p=tmp_p*d(j)/2.;
                w=cross(sv,ref_vec_next);w=w/norm(w);
                [tmp_pointnumber,tmp]=size(AllPoint);
                
                tmp_label=ones(m1,1);
                if index_sec_pt==2
                    tmp_label=tip_label*tmp_label;
                    tip_label=tip_label+1;
                else
                    tmp_label=in_point_label*tmp_label;
                end
                tmp_label(boundary_point_index,:)=wall_label;
                
                for inode=1:m1
                    tmp_p(inode,:)= tmp_p(inode,:)+location(j,:)-location(1,:);
                end
                for ii=1:m1
                    tmp_veloctity(ii,1:3)=sv/norm(sv)*norm(velocity_value(ii));
                end
                AllVelocity=[AllVelocity; tmp_veloctity];
                AllPoint=[AllPoint;tmp_p];
                AllLabel=[AllLabel;tmp_label];
                for k=1:m2
                    tmp_e(k,1:4)=template_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
                end
                LayerElement=[LayerElement;tmp_e(:,1:4)];
                NodeLayer{j,1}=tmp_p;
                NodeLayer{j,2}=tmp_e(:,1:4);
                NodeLayer{j,3}=tmp_e(:,1:4);
            end
            continue;
        end
        if(termination(index_end_pt))
            n_insert_layer=n_insert_layer+1;
            sec_vec_start=Segment_Vector{index_start_pt,id(index_start_pt)};
            ref_vec_start=NodeLayer{index_start_pt,4};
            w=cross(sec_vec_start,ref_vec_start); w=w/norm(w);
            ref_vec_next=ref_vec_start/norm(ref_vec_start);
            
            for index_sec_pt=2:sect_ptnum
                j=bif_term_pt{index_sec}(index_sec_pt-1);
                i=bif_term_pt{index_sec}(index_sec_pt);
                sv=Segment_Vector{i,j};
                ref_vec_next=cross(w,sv); ref_vec_next=ref_vec_next/norm(ref_vec_next);
                w=cross(sv,ref_vec_next);w=w/norm(w);
                for ii=1:m1
                    tmp_p(ii,:)=template_p(ii,1)*ref_vec_next+template_p(ii,2)*w;
                end
                tmp_p=tmp_p*d(i)/2.;
                
                [tmp_pointnumber,tmp]=size(AllPoint);
                
                tmp_label=ones(m1,1);
                if index_sec_pt==sect_ptnum
                    tmp_label=tip_label*tmp_label;
                    tip_label=tip_label+1;
                else
                    tmp_label=in_point_label*tmp_label;
                end
                tmp_label(boundary_point_index,:)=wall_label;
                
                for inode=1:m1
                    tmp_p(inode,:)= tmp_p(inode,:)+location(i,:)-location(1,:);
                end
                for ii=1:m1
                    tmp_veloctity(ii,1:3)=sv/norm(sv)*norm(velocity_value(ii));
                end
                AllVelocity=[AllVelocity; tmp_veloctity];
                AllPoint=[AllPoint;tmp_p];
                AllLabel=[AllLabel;tmp_label];
                for k=1:m2
                    tmp_e(k,1:4)=template_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
                end
                LayerElement=[LayerElement;tmp_e(:,1:4)];
                NodeLayer{i,1}=tmp_p;
                NodeLayer{i,2}=tmp_e(:,1:4);
                NodeLayer{i,3}=tmp_e(:,1:4);
            end
            continue;
        end
        % bif-bif segments
        sec_vec_start=Segment_Vector{index_start_pt,id(index_start_pt)};
        sec_vec_end=Segment_Vector{index_end_pt,id(index_end_pt)};
        ref_vec_start=NodeLayer{index_start_pt,4};
        ref_vec_end=NodeLayer{index_end_pt,4};
        
        rotate_axis=cross(ref_vec_start,ref_vec_end);
        angle_total=acos(dot(ref_vec_start,ref_vec_end)/norm(ref_vec_start)/norm(ref_vec_end));
        if angle_total>pi/4 && angle_total<=3*pi/4
            angle_total=angle_total-pi/2;
            if dot(rotate_axis,sec_vec_end)>0
                for index_ele=1:m2/4
                    tmp_end_element(index_ele+m2/4*0,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*3,[3 4 1 2]);
                    tmp_end_element(index_ele+m2/4*1,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*0,:);
                    tmp_end_element(index_ele+m2/4*2,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*1,[3 4 1 2]);
                    tmp_end_element(index_ele+m2/4*3,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*2,:);
                end
                NodeLayer{index_end_pt,3}=tmp_end_element;
            elseif dot(rotate_axis,sec_vec_end)<0
                for index_ele=1:m2/4
                    tmp_end_element(index_ele+m2/4*0,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*1,:);
                    tmp_end_element(index_ele+m2/4*1,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*2,[3 4 1 2]);
                    tmp_end_element(index_ele+m2/4*2,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*3,:);
                    tmp_end_element(index_ele+m2/4*3,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*0,[3 4 1 2]);
                end
                NodeLayer{index_end_pt,3}=tmp_end_element;
            end
        elseif angle_total>3*pi/4 && angle_total<=pi
            angle_total=angle_total-pi;
            for index_ele=1:m2/4
                tmp_end_element(index_ele+m2/4*0,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*2,[3 4 1 2]);
                tmp_end_element(index_ele+m2/4*1,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*3,[3 4 1 2]);
                tmp_end_element(index_ele+m2/4*2,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*0,[3 4 1 2]);
                tmp_end_element(index_ele+m2/4*3,:)=NodeLayer{index_end_pt,2}(index_ele+m2/4*1,[3 4 1 2]);
            end
            NodeLayer{index_end_pt,3}= tmp_end_element;
        end
        
        angle_per=angle_total/(n_insert_layer+1);
        
        fprintf('sec_index=%d angle_total=%f angle_per=%f\n',index_sec,angle_total*180/pi,angle_per*180/pi);
        
        w=cross(sec_vec_start,ref_vec_start); w=w/norm(w);
        ref_vec_next=ref_vec_start/norm(ref_vec_start);
        for index_sec_pt=2:sect_ptnum-1
            j=bif_term_pt{index_sec}(index_sec_pt-1);
            i=bif_term_pt{index_sec}(index_sec_pt);
            sv=Segment_Vector{i,j};
            
            ref_vec_next=cross(w,sv); ref_vec_next=ref_vec_next/norm(ref_vec_next);
            w=cross(sv,ref_vec_next);w=w/norm(w);
            for ii=1:m1
                tmp_p(ii,:)=template_p(ii,1)*ref_vec_next+template_p(ii,2)*w;
            end
            tmp_p=tmp_p*d(i)/2.;
            
            
            if angle_total~=0
                tmp_p=RotateAroundAxis(tmp_p,sv,angle_per*(index_sec_pt-1));
            else
                tmp_p=tmp_p;
            end
            
            
            for ii=1:m1
                tmp_p(ii,:)=template_p(ii,1)*ref_vec_next+template_p(ii,2)*w;
            end
            tmp_p=tmp_p*d(i)/2.;
            
            [tmp_pointnumber,tmp]=size(AllPoint);
            
            tmp_label=ones(m1,1);
            tmp_label=in_point_label*tmp_label;
            tmp_label(boundary_point_index,:)=wall_label;
            
            for inode=1:m1
                tmp_p(inode,:)= tmp_p(inode,:)+location(i,:)-location(1,:);
            end
            for ii=1:m1
                tmp_veloctity(ii,1:3)=sv/norm(sv)*norm(velocity_value(ii));
            end
            AllVelocity=[AllVelocity; tmp_veloctity];
            AllPoint=[AllPoint;tmp_p];
            AllLabel=[AllLabel;tmp_label];
            for k=1:m2
                tmp_e(k,1:4)=template_e(k,2:5)+tmp_pointnumber*[1,1,1,1];
            end
            LayerElement=[LayerElement;tmp_e(:,1:4)];
            NodeLayer{i,1}=tmp_p;
            NodeLayer{i,2}=tmp_e(:,1:4);
            NodeLayer{i,3}=tmp_e(:,1:4);
        end
    end
else
end

%% Calculate the points and index of the layers in each segment
for index_sec=1:n_sect
    [sect_ptnum,tmp]=size(sect_pt{index_sec}');
    for index_sec_pt=1:sect_ptnum-1
        j=sect_pt{index_sec}(index_sec_pt);
        i=sect_pt{index_sec}(index_sec_pt+1);
        segment_layer=n_layer(i,j);
        sv=Segment_Vector{i,j};
        if(branch(j))
            sv=Segment_Vector{i,j};
            if(ijk_label(i)==1)
                tmp_start_e=NodeLayer{j,3};
                tmp_start_p=BranchLayerPoint{j,2};
            elseif(ijk_label(i)==2)
                tmp_start_e=NodeLayer{j,4};
                tmp_start_p=BranchLayerPoint{j,3};
            end
        else
            tmp_start_p=NodeLayer{j,1};
            tmp_start_e=NodeLayer{j,2};
        end
        if(branch(i))
            tmp_end_p=BranchLayerPoint{i,1};
            tmp_end_e=NodeLayer{i,2};
        else
            tmp_end_p=NodeLayer{i,1};
            tmp_end_e=NodeLayer{i,3};
        end
        if segment_layer==1
            [tmp_pointnumber,tmp]=size(AllPoint);
            AllElement=[AllElement;tmp_start_e,tmp_end_e];
        else
            for k=1:segment_layer-1
                [tmp_pointnumber,tmp]=size(AllPoint);
                [element_index,tmp]=size(AllElement);
                [layer_index,tmp]=size(LayerElement);
                
                % record the element information around bifurcation
                if(branch(i) && k==segment_layer-1)
                    index_bif=find(bif_pt(:,1)==i);
                    bif_ele(index_bif,2)=tmp_pointnumber;
                end
                if(branch(j) && k==1)
                    index_bif=find(bif_pt(:,1)==j);
                    if(ijk_label(i)==1 )
                        bif_ele(index_bif,3)=tmp_pointnumber;
                    elseif(ijk_label(i)==2)
                        bif_ele(index_bif,4)=tmp_pointnumber;
                    end
                end
                
                for inode=1:m1
                    tmp_p(inode,:)=(tmp_end_p(inode,:)*k+tmp_start_p(inode,:)*(segment_layer-k))/segment_layer;
                end
                
                tmp_label=ones(m1,1);
                tmp_label=in_point_label*tmp_label;
                tmp_label(boundary_point_index,:)=wall_label;
                
                for ii=1:m1
                    tmp_veloctity(ii,1:3)=sv/norm(sv)*norm(velocity_value(ii));
                end
                AllVelocity=[AllVelocity; tmp_veloctity];
                
                AllPoint=[AllPoint;tmp_p];
                AllLabel=[AllLabel;tmp_label];
                for iele=1:m2
                    tmp_e(iele,1:4)=template_e(iele,2:5)+tmp_pointnumber*[1,1,1,1];
                end
                LayerElement=[LayerElement;tmp_e(:,1:4)];
                if(k==1)
                    AllElement=[AllElement;tmp_start_e,tmp_e(:,1:4)];
                    AllElement=[AllElement;tmp_e(:,1:4),zeros(m2,4)];
                else
                    AllElement((element_index+1-m2):element_index,5:8)=tmp_e(:,1:4);
                    AllElement=[AllElement;tmp_e(:,1:4),zeros(m2,4)];
                end
            end
            [element_index,tmp]=size(AllElement);
            AllElement((element_index+1-m2):element_index,5:8)=tmp_end_e;
        end
    end
end

% Deal with bifurcation point
for index_bif=1:n_bif
    startpt_b=bif_ele(index_bif,1);
    startpt_i=bif_ele(index_bif,2);
    startpt_j=bif_ele(index_bif,3);
    startpt_k=bif_ele(index_bif,4);
    
    [num_bc,tmp]=size(leftbc_index);
    for i=1:num_bc
        AllPoint(startpt_b+botbc_index(i),:)=PointAlign(AllPoint(startpt_j+rightbc_index(i),:),AllPoint(startpt_k+leftbc_index(i),:),AllPoint(startpt_b+botbc_index(i),:));
        AllPoint(startpt_b+leftbc_index(i),:)=PointAlign(AllPoint(startpt_i+leftbc_index(i),:),AllPoint(startpt_j+leftbc_index(i),:),AllPoint(startpt_b+leftbc_index(i),:));
        AllPoint(startpt_b+rightbc_index(i),:)=PointAlign(AllPoint(startpt_i+rightbc_index(i),:),AllPoint(startpt_k+rightbc_index(i),:),AllPoint(startpt_b+rightbc_index(i),:));
    end
    for i=1:2
        AllPoint(startpt_b+extra_pt(i,1),:)=Projection(AllPoint(startpt_b+extra_pt(i,2),:),AllPoint(startpt_b+extra_pt(i,3),:),AllPoint(startpt_b+extra_pt(i,4),:),AllPoint(startpt_b+extra_pt(i,1),:));
        AllPoint(startpt_i+extra_pt(i,1),:)=Projection(AllPoint(startpt_b+extra_pt(i,2),:),AllPoint(startpt_b+extra_pt(i,3),:),AllPoint(startpt_b+extra_pt(i,4),:),AllPoint(startpt_i+extra_pt(i,1),:));
        AllPoint(startpt_j+extra_pt(i,1),:)=Projection(AllPoint(startpt_b+extra_pt(i,2),:),AllPoint(startpt_b+extra_pt(i,3),:),AllPoint(startpt_b+extra_pt(i,4),:),AllPoint(startpt_j+extra_pt(i,1),:));
        AllPoint(startpt_k+extra_pt(i,1),:)=Projection(AllPoint(startpt_b+extra_pt(i,2),:),AllPoint(startpt_b+extra_pt(i,3),:),AllPoint(startpt_b+extra_pt(i,4),:),AllPoint(startpt_k+extra_pt(i,1),:));
    end
end
%% Output Mesh
%write the vtk file for view
[n_point,ff]=size(AllPoint);
[n_element,ff]=size(AllElement);
[n_element_layer,ff]=size(LayerElement);
fid1=fopen(velocity_output,'w');
for ii=1:n_point
    fprintf(fid1,'%f %f %f\n',AllVelocity(ii,1:3));
end

%write the mesh vtk file
fid2=fopen(hex_output,'w');
fprintf(fid2,'%s\n','# vtk DataFile Version 3.1 ');
fprintf(fid2,'%s\n','for LSEConsole');
fprintf(fid2,'%s\n','ASCII');
fprintf(fid2,'%s\n','DATASET UNSTRUCTURED_GRID');
fprintf(fid2,'%s %d %s\n','POINTS',n_point,'FLOAT');
for ii=1:n_point
    fprintf(fid2,'%f %f %f\n',AllPoint(ii,1:3));
end
fprintf(fid2,'%s %d %d\n','CELLS',n_element,9*n_element);
for  ii=1:n_element
    fprintf(fid2,'%d %d %d %d %d %d %d %d %d\n', 8 ,AllElement(ii,:));
end
fprintf(fid2,'%s %d\n','CELL_TYPES',n_element);
for  ii=1:n_element
    fprintf(fid2,'%d\n',12);
end
fprintf(fid2,'%s %d\n','POINT_DATA',n_point);
fprintf(fid2,'%s\n','SCALARS label float 1');
fprintf(fid2,'%s\n','LOOKUP_TABLE default');
for  ii=1:n_point
    fprintf(fid2,'%d\n',AllLabel(ii,:));
end
fclose('all');