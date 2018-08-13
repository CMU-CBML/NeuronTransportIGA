clc;
clear;
addpath(genpath(pwd));
start_trees;
%% User settings
% Set smooth parameters
n_noisesmooth=150;      % set iteration steps for noise smooth
ratio_bifur_node=0.01;  % set bifurcation nodes smooth ratio
ratio_noisesmooth=0.01; % set noise smooth ratio
seg_length=0.5;         % set bezier smooth segments length 

% input and output path setting
io_path='..//io//bifurcation1//'; % user set the input and output path 
neuron_name='bifurcation1';       % user set the input skeleton file name 
    %neuron_name='cell3traceRN1';
    %neuron_name='nelson2';
    %neuron_name='purkinje_modify3';

input_file=[io_path,neuron_name,'.swc'];
smooth_file=[io_path,neuron_name,'_smooth.swc'];
tangent_file=[io_path,neuron_name,'_tangent.txt'];

trees{1}=load_tree(input_file);
%% Extract coordinates and connection information
tangent_vec=[];
location=[trees{1}.X,trees{1}.Y,trees{1}.Z];
d=trees{1}.D;
id=idpar_tree';
id(1)=0;
ipar=ipar_tree();

branch_input=B_tree(trees{1});
termination=T_tree(trees{1});
[sect, vec]=dissect_tree(trees{1});
[n_sect, tmp]= size(sect);
[n_bif,tmp]=size(find(branch_input==1));
sect_point=cell(n_sect,1);
bif_term_pt=cell(n_sect,1);
bif_pt=cell(n_bif,1);

[nx,ny]=size(id);
[m,n]=size(find(trees{1}.dA~=0));
[mm,nn]=size(ipar);

OutputData=[];

%% Extract all points for each branch and all bifurcation points
for i=1:n_sect
    tmp_par=ipar(sect(i,2),:);
    loc=find(tmp_par==sect(i,1));
    sect_point{i}=tmp_par(loc:-1:1);
end

index_bif=1;
for sec_index=1:n_sect
    [sect_ptnum,tmp]=size(sect_point{sec_index}');
    start_pt_index=sect_point{sec_index}(1);
    end_pt_index=sect_point{sec_index}(sect_ptnum);
    
    bif_term_start=1;
    bif_term_end=sect_ptnum;
    
    if(branch_input(end_pt_index))
        par_node=id(end_pt_index);
        ans_tmp=find(trees{1}.dA(:,end_pt_index)==1);
        bifur1_node=ans_tmp(1);
        bifur2_node=ans_tmp(2);
        bif_pt{index_bif}=[end_pt_index,par_node,bifur1_node,bifur2_node];
        index_bif=index_bif+1;
        bif_term_end=sect_ptnum-1;
    end
end


%% Noise Smooth
for index_smooth=1:n_noisesmooth
% Optimize bifurcation point to avoid bad geometry during mesh generation
% Move the bifurcation points to the middle of three neighbour nodes
    for ii=1:index_bif-1
        pt_b=[trees{1}.X(bif_pt{ii}(1)),trees{1}.Y(bif_pt{ii}(1)),trees{1}.Z(bif_pt{ii}(1))];  d_b=trees{1}.D(bif_pt{ii}(1));
        pt_i=[trees{1}.X(bif_pt{ii}(2)),trees{1}.Y(bif_pt{ii}(2)),trees{1}.Z(bif_pt{ii}(2))];  d_i=trees{1}.D(bif_pt{ii}(2));
        pt_j=[trees{1}.X(bif_pt{ii}(3)),trees{1}.Y(bif_pt{ii}(3)),trees{1}.Z(bif_pt{ii}(3))];  d_j=trees{1}.D(bif_pt{ii}(3));
        pt_k=[trees{1}.X(bif_pt{ii}(4)),trees{1}.Y(bif_pt{ii}(4)),trees{1}.Z(bif_pt{ii}(4))];  d_k=trees{1}.D(bif_pt{ii}(4));
        pt_b_after=pt_b*(1-ratio_bifur_node)+(pt_i+pt_j+pt_k)/3*ratio_bifur_node;
        d_b_after=d_b*(1-ratio_bifur_node)+(d_i+d_j+d_k)/3*ratio_bifur_node;
        trees{1}.X(bif_pt{ii}(1))=pt_b_after(1);
        trees{1}.Y(bif_pt{ii}(1))=pt_b_after(2);
        trees{1}.Z(bif_pt{ii}(1))=pt_b_after(3);
        trees{1}.D(bif_pt{ii}(1))=d_b_after;
    end
% Eliminate noise point
    for i=1:n_sect
        AA=[trees{1}.X(sect_point{i}),trees{1}.Y(sect_point{i}),trees{1}.Z(sect_point{i}), trees{1}.D(sect_point{i})];
        BB=NoiseSmooth(AA,ratio_noisesmooth,1);
        trees{1}.X(sect_point{i})=BB(:,1);
        trees{1}.Y(sect_point{i})=BB(:,2);
        trees{1}.Z(sect_point{i})=BB(:,3);
        trees{1}.D(sect_point{i})=BB(:,4);
    end
end

location=[trees{1}.X,trees{1}.Y,trees{1}.Z];
d=trees{1}.D;

inter_pt=[];
for i=1:n_sect
    inter_pt=[inter_pt,sect_point{i}(2:end-1)];
end
trees{2}=delete_tree(trees{1},inter_pt);
branch2=B_tree(trees{2});
termination2=T_tree(trees{2});

[sect_after, vec_after]=dissect_tree(trees{2});
location_after=[trees{2}.X,trees{2}.Y,trees{2}.Z];
d_after=trees{2}.D;
[ptnum,tmp]=size(d_after);
branch_insert_vector=cell(ptnum,3);
%% Bezier Smooth
for i=1:n_sect
    start_index=sect_after(i,1);
    end_index=sect_after(i,2);
    if(branch2(end_index) && branch2(start_index))
        mode=1;
    elseif(termination2(end_index))
        mode=2;
    elseif(start_index==1)
        mode=3;
    end
    [tmp_XYZ,tmp_D, tmp_tangent]=BsplineSmooth(location(sect_point{i},:),d(sect_point{i}),seg_length ,mode);
    [n_insert,tmp]=size(tmp_D);
    tangent_vec=[tangent_vec;start_index tmp_tangent(1,:)];
    for j=2:n_insert-1
        [index,tmp]=size(trees{2}.R);
        index=index+1;
        tangent_vec=[tangent_vec;index tmp_tangent(j,:)];
        if (j==2)
            trees{2}=insert_tree(trees{2},[index,2,tmp_XYZ(j,:), tmp_D(j),sect_after(i,1)]);
        else
            trees{2}=insert_tree(trees{2},[index,2,tmp_XYZ(j,:), tmp_D(j),index-1]);
        end
    end
    trees{2}=recon_tree(trees{2},sect_after(i,2),index,'none');
    tangent_vec=[tangent_vec;end_index tmp_tangent(end,:)];
end


[tmp_C, ia, ic]=unique(tangent_vec(:,1),'sorted');
output_tangent=tangent_vec(ia,:);

%% Output and Visualization
trees{1}.R(:,1)=2;
trees{2}.R(:,1)=2;

swc_tree(trees{2},smooth_file);
fid4=fopen(tangent_file,'w');
[n_vec,tmp]=size(output_tangent);
for ii=1:n_vec
    fprintf(fid4,'%f %f %f\n',output_tangent(ii,2:4));
end
fclose('all');

figure(1);clf;xplore_tree(trees{1})
figure(2);clf;xplore_tree(trees{2})