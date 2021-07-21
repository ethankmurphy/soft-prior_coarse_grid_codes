%-------------------------------------------------------------------------------
%
% This function constructs coarse subgrids based on the underlying FEM
% mesh, the interior surfaces, and flags indicating any outer or inner
% surfaces. All the surfaces are setup so that the normals are pointing
% outward.
%
%
% msh - is our standard mesh structure
%   msh.node - nodes
%   msh.elem - 3D tetrahedra
%   msh.face - boundary triangles
%   msh.elec - list of electrodes
% intfuncs   - cell array of interior functions (could be 1 or multiple
%              interior regions (an interior point is negative)
% outfunc    - function describing the outer boundary of the sub-region
%              (an interior point is negative)
% dx         - size of coarse grid node spacing
% nscl       - scaling of nodes
% dbg_flg    - flag
%
%-------------------------------------------------------------------------------
function [cgrid,c2f] = constr_coarsesub_grid_genfuncs(msh,intfuncs,outfunc,dx,nscl,dbg_flg,cgrid)                                                      
%-------------------------------------------------------------------------------
% nscl    = 1000;
nds     = msh.node*nscl;

%-------------------------------------------------------------------------------
% Find FEM mesh nodes that are within the domain
% Find the interior nodes - these nodes will be removed
is_int = [];
for n = 1:length(intfuncs)
    is_int = [is_int; find( intfuncs{n}(nds) < 0)];
end
% Find the subregion nodes
if isempty(outfunc)
    isr = setdiff( (1:size(nds,1))',is_int);
else
    isr = find( outfunc(nds) < 0);
    isr = setdiff(isr,is_int);
end

%-------------------------------------------------------------------------------
if dbg_flg == 1
    figure
    hold on
    istmp1 = find( nds(isr,2) > 0);    
    istmp2 = find( nds(is_int,1) > 0);    
    plot3(nds(isr(istmp1),1),nds(isr(istmp1),2),nds(isr(istmp1),3),'.k','markersize',4)
    plot3(nds(is_int(istmp2),1),nds(is_int(istmp2),2),nds(is_int(istmp2),3),'.r','markersize',8)
    legend('Outer','Inner')
    title('Checking FEM nodes in and outside of the domain')
    lbl_fmt_fig('X','Y','','','Z',12)
end

%-------------------------------------------------------------------------------
% Construct a coarse grid that lies within the subdomain
if nargin < 7
    cgrid = construct_standard_rectgrid(isr,nds,msh.elem,dx);
    
    if dbg_flg == 1
        disp('Finished: Standard Rectangular grid')
    end
    %-------------------------------------------------------------------------------
    % Restrict coarse grid to the given subdomain
    is_int_c = [];
    for n = 1:length(intfuncs)
        is_int_c = [is_int_c; find( intfuncs{n}(cgrid.node) < 0)];
    end
    % Find the subregion nodes
    if isempty(outfunc)
        isc = setdiff( (1:size(cgrid.node,1))',is_int);
    else
        isc = find( outfunc(cgrid.node) < 0);
        isc = setdiff(isc,is_int_c);
    end
    if dbg_flg == 1
        disp('Finished: Coarse grid subregion')
        cgrid
        %---------------------------------------------------------------------------
        figure
        hold on
        plot3(cgrid.node(isc,1),cgrid.node(isc,2),cgrid.node(isc,3),'.k','markersize',4)
        plot3(cgrid.node(is_int_c,1),cgrid.node(is_int_c,2),cgrid.node(is_int_c,3),'.r','markersize',8)
        legend('Outer','Inner')
        title('Checking Coarse nodes in and outside of the domain')
        lbl_fmt_fig('X','Y','','','Z',12)
    end
    
    %-------------------------------------------------------------------------------
    % Keep only subdomain points and scale to meters
    cgrid.node = cgrid.node(isc,:)/nscl; % The rest of the calc assumes in meters for everything
    cgrid.dx   = cgrid.dx/nscl; % The rest of the calc assumes in meters for everything
    cgrid.dy   = cgrid.dy/nscl; % The rest of the calc assumes in meters for everything
    cgrid.dz   = cgrid.dz/nscl; % The rest of the calc assumes in meters for everything
end
cgrid.isr  = isr;
%-------------------------------------------------------------------------------
% Construct the coarse to fine mapping (fine is the full FEM grid)
num_nod_f = size(msh.node,1);
num_nod_c = size(cgrid.node,1);
c2f       = spalloc(num_nod_f,num_nod_c,num_nod_f);

%-------------------------------------------------------------------------------
% 2.  Cycle through the coarse nodes and find the closest fine (FEM) node that
%     is within the given domain
for i=1:num_nod_c
    
    %---------------------------------------------------------------------------
    % squared distance of the i-th fine node to any interior coarse node
    distance=( ...
        (msh.node(isr,1)-cgrid.node(i,1)).^2 + ...
        (msh.node(isr,2)-cgrid.node(i,2)).^2 + ...
        (msh.node(isr,3)-cgrid.node(i,3)).^2);
    
    %---------------------------------------------------------------------------
    % find index to closest coarse node
    [md,ind]        = min(distance); % find the closest corse node
    if md < 1e-10
        md = 1e-10;
        disp('got in here')
    end
    c2f(isr(ind),i) = 1/md;
    
end % for i
if dbg_flg == 1
    disp('Finished: c2f over coarse grid nodes')
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% We want average the the fine FEM node value by all the assigned coarse
% nodes relative to their distance from the fine FEM node
% c2f is size (num. fine) x (num. coarse)
num_assigned_f   = sum(double(c2f(isr,:)>1e-12),2); % Number of assigned coarse
% nodes to each fine node
%-------------------------------------------------------------------------------
% For those FEM nodes with only 1 assigned coarse node. Set the value equal
% to 1.
is                  = find(num_assigned_f == 1);
[max_percol,imxs]   = max( c2f(isr,:),[],2);
for n = 1:length(is)
    c2f( isr(is(n)),: ) = 0;
    c2f( isr(is(n)),imxs(is(n)) ) = 1;
end
fnd_lb     = zeros(num_nod_f,1);
fnd_lb(isr(is)) = 1;
%-------------------------------------------------------------------------------
% Loop through the remaining fine nodes with multiple coarse nodes
% assigned. Construct a weighted average approximation
is                  = find(num_assigned_f > 1);
sum_percol          = sum( c2f(isr,:),2);
for n = 1:length(is)
    c2f( isr(is(n)),: ) = 1/sum_percol(is(n))*c2f( isr(is(n)),: );
end
if dbg_flg == 1
    disp('Finished: c2f fixes')
end
fnd_lb(isr(is)) = 2;
%-------------------------------------------------------------------------------
% Loop through the unassociated fine nodes
is                  = find(num_assigned_f == 0);
if dbg_flg == 1
    disp(['Num. unassociated fine nodes: ',num2str(length(is))])
end
fnd_lb(isr(is)) = 3;
for i = 1:length(is)
    %---------------------------------------------------------------------------
    % squared distance of the i-th fine node to any interior coarse node
    distance=( ...
        (msh.node(isr(is(i)),1)-cgrid.node(:,1)).^2 + ...
        (msh.node(isr(is(i)),2)-cgrid.node(:,2)).^2 + ...
        (msh.node(isr(is(i)),3)-cgrid.node(:,3)).^2);
    
    %---------------------------------------------------------------------------
    % find index to closest coarse node
    [md,ind]=min(distance); % find the closest corse node
    c2f(isr(is(i)),ind)= 1;
    
end % for i
if dbg_flg == 1
    disp('Finished: c2f unassociated fine nodes')
end
cgrid.fnd_lb = fnd_lb;

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%
%                       get_nds_inds_within_subreg
%
% This function gets the indices of the nodes (that come wherever) that are
% within the given subregion defined by inner and outer surfaces.
%
%-------------------------------------------------------------------------------
function cgrid = construct_standard_rectgrid(isr,nds,elem,dx)
whos
ris = setdiff((1:size(nds,1))',unique(isr));
[nds,elem] = update_tris_due_to_rmvnodes(nds,elem,ris);
whos
% is_crazy1 = find( elem(:,1) == 0);
% if ~isempty(is_crazy1)
%     elem(is_crazy1,:)
% end


%-------------------------------------------------------------------------------
% Get the boundary of the coarse grid
xmax = max(nds(:,1)); % finds x max of finer msh
xmin = min(nds(:,1));
ymax = max(nds(:,2)); % finds y max of finer msh
ymin = min(nds(:,2));

zmax = max(nds(:,3)); % finds y max of finer msh
zmin = min(nds(:,3));


% Set dx, dy and dz, and compute node positions
if length(dx) == 1
    dx = repmat(dx,1,3);
end

%-------------------------------------------------------------------------------
% X:
numpix(1)= max([round((xmax-xmin)/dx(1)) 3]);
x        = linspace(xmin+dx(1)/2,xmax-dx(1)/2,numpix(1));
cgrid.dx = x(2)-x(1);
%-------------------------------------------------------------------------------
% X:
numpix(2)= max([round((ymax-ymin)/dx(2)) 3]);
y        = linspace(ymin+dx(2)/2,ymax-dx(2)/2,numpix(2));
cgrid.dy = y(2)-y(1);
%-------------------------------------------------------------------------------
% X:
numpix(3)= max([round((zmax-zmin)/dx(3)) 3]);
z        = linspace(zmin+dx(3)/2,zmax-dx(3)/2,numpix(3));
cgrid.dz = z(2)-z(1);

%-------------------------------------------------------------------------------
% generate the corse grid
[X,Y,Z]=meshgrid(x,y,z);

% rearrange the grid into a (n x 3) matrix of nodes
cgrid.node = [reshape(X,numpix(1)*numpix(2)*numpix(3),1) reshape(Y,numpix(1)*numpix(2)*numpix(3),1) reshape(Z,numpix(1)*numpix(2)*numpix(3),1)];

% lets see which corse grid nodes are falling within the fine mesh
[tind,p]=tsearchn(nds, elem, cgrid.node); % we will have NaNs for nodes outside the mesh
ind=find(~isnan(p(:,1))); % find an idex to interior corse nodes (not equal to NaN) - we look into p instead of tind,
% as tind returns NaNs for tetra outside the convex hull, while p, the matrix of interpolation weights, returns NaNs for points outside the mesh

% lets keep only the interior nodes
cgrid.node=cgrid.node(ind,:);

