%--------------------------------------------------------------------------
%
% This codes gives two examples of
%
%--------------------------------------------------------------------------
clear
clc
close all

%--------------------------------------------------------------------------
test_typ = 1;

%--------------------------------------------------------------------------
% Load the mesh with encoded surfaces
if test_typ == 1
    mshfnme = 'sphere_1inc_r10_2_h0o8_0o3_0o6';
elseif test_typ == 2
    mshfnme = 'sphere_shell_1inc_r10_8_2_h0o8_0o3_0o6_0o6';
end
eval(['load dat/',mshfnme,' msh out_p out_f in_ps in_fs eps ediam'])
msh.node = msh.node*1000;

%--------------------------------------------------------------------------
% Plot the mesh with its interior surfaces
figure;hold on
simpplot(msh.node,msh.elem,'p(:,2)>0');hold on
pcols = {'b','m','c'};
trisurf(out_f,out_p(:,1),out_p(:,2),out_p(:,3),'linestyle','none','facecolor',pcols{1},'facealpha',0.3)
for n = 1:length(in_ps)
    trisurf(in_fs{n},in_ps{n}(:,1),in_ps{n}(:,2),in_ps{n}(:,3),'linestyle','none','facecolor',pcols{n+1},'facealpha',0.3)
end
plot_msh_elecs_only(msh,1,0,1)
axis equal
axis(10*[-1 1 -1 1 -1 1])
camlight left
view(3)
title('FEM mesh with subdomains')

%--------------------------------------------------------------------------
if test_typ == 1
    %----------------------------------------------------------------------
    % One inclusion case
    %----------------------------------------------------------------------
    % Set the functions describing each domain/subdomain. The function 
    % should be negative inside the domain and positve outside. The 
    % surface triangulations in the mesh mat-file are setup to have 
    % outward pointing normals.
    incfunc{1} = @(p) (eval_signed_dist(in_ps{1},in_fs{1},[],p,0));
     
    dx_out = 2; % Voxel size in outer portion
    dx_in  = 1; % Voxel size in inner portion
    %----------------------------------------------------------------------
    % Construct the coarse grid for the outer portion of the domain
    [cgrid_out,c2f_out] = constr_coarsesub_grid_genfuncs(msh,incfunc,[],dx_out,1,0);
    % Make a standard regularization matrix in the subdomain
    if size(cgrid_out.node,1) == 1
        Lout = 1;
    else
        Lout = ndrm_laplacian_regmatrix(cgrid_out);
    end
    %----------------------------------------------------------------------
    % Construct the coarse grid for the outer portion of the domain
    [cgrid_in,c2f_in] = constr_coarsesub_grid_genfuncs(msh,[],incfunc{1},dx_in,1,0);
    % Make a standard regularization matrix in the subdomain
    if size(cgrid_in.node,1) == 1 
        Lin = 1;
    else
        Lin = ndrm_laplacian_regmatrix(cgrid_in);
    end
    %----------------------------------------------------------------------
    figure;hold on
    plot3(cgrid_out.node(:,1),cgrid_out.node(:,2),cgrid_out.node(:,3),'.k','markersize',12)
    plot3(cgrid_in.node(:,1),cgrid_in.node(:,2),cgrid_in.node(:,3),'.r','markersize',16)
    title('Each coarse grid')
    %----------------------------------------------------------------------
    % Put the coarse grids, coarse to fine maps, and regularization 
    % matrices into cell arrays to a merging function can be used
    cgrids = {cgrid_out,cgrid_in};
    c2fs   = {c2f_out,c2f_in};
    Ls     = {Lout,Lin};
    
elseif test_typ == 2
    %----------------------------------------------------------------------
    % One shell and an inclusion in the inner portion
    %----------------------------------------------------------------------
    % Set the functions describing each domain/subdomain. The function 
    % should be negative inside the domain and positve outside. The 
    % surface triangulations in the mesh mat-file are setup to have 
    % outward pointing normals.
    incfunc{1}   = @(p) (eval_signed_dist(in_ps{2},in_fs{2},[],p,0));
    brainfunc{1} = @(p) (eval_signed_dist(in_ps{1},in_fs{1},[],p,0));
     
    dx_shlout = 2; % Voxel size in outer, shell portion (most outer)
    dx_out    = 1.5; % Voxel size in outer portion
    dx_in     = 1; % Voxel size in inner portion
    %----------------------------------------------------------------------
    % Construct the coarse grid for the outer-most portion of the domain
    [cgrid_out1,c2f_out1] = constr_coarsesub_grid_genfuncs(msh,brainfunc,[],dx_shlout,1,0);
    % Make a standard regularization matrix in the subdomain
    if size(cgrid_out1.node,1) == 1
        Lout1 = 1;
    else
        Lout1 = ndrm_laplacian_regmatrix(cgrid_out1);
    end
    %----------------------------------------------------------------------
    % Construct the coarse grid for the inner large portion (maybe brain) 
    [cgrid_out2,c2f_out2] = constr_coarsesub_grid_genfuncs(msh,incfunc,brainfunc{1},dx_shlout,1,0);
    % Make a standard regularization matrix in the subdomain
    if size(cgrid_out2.node,1) == 1
        Lout2 = 1;
    else
        Lout2 = ndrm_laplacian_regmatrix(cgrid_out2);
    end
    %----------------------------------------------------------------------
    % Construct the coarse grid for the outer portion of the domain
    [cgrid_in,c2f_in] = constr_coarsesub_grid_genfuncs(msh,[],incfunc{1},dx_in,1,0);
    % Make a standard regularization matrix in the subdomain
    if size(cgrid_in.node,1) == 1
        Lin = 1;
    else
        Lin = ndrm_laplacian_regmatrix(cgrid_in);
    end
    
    %----------------------------------------------------------------------
    figure;hold on
    plot3(cgrid_out1.node(:,1),cgrid_out1.node(:,2),cgrid_out1.node(:,3),'.k','markersize',12)
    plot3(cgrid_out2.node(:,1),cgrid_out2.node(:,2),cgrid_out2.node(:,3),'.b','markersize',12)
    plot3(cgrid_in.node(:,1),cgrid_in.node(:,2),cgrid_in.node(:,3),'.r','markersize',16)
    view(3)
    title('Each coarse grid')
    %----------------------------------------------------------------------
    % Put the coarse grids, coarse to fine maps, and regularization 
    % matrices into cell arrays to a merging function can be used
    cgrids = {cgrid_out1,cgrid_out2,cgrid_in};
    c2fs   = {c2f_out1,c2f_out2,c2f_in};
    Ls     = {Lout1,Lout2,Lin};
    
end

%--------------------------------------------------------------------------
% Merge the coarse grids, coarse to fine maps, and regularization matrices
[cnds,c2f,L] = merge_grids_maps_Lmats(cgrids,c2fs,Ls);

%--------------------------------------------------------------------------
% Check everything is working... 
check_softprior_cgrid_c2f_L(cnds,c2f,L,msh);
