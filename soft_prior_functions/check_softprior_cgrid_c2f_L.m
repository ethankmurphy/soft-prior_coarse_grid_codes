%--------------------------------------------------------------------------
% 
% Check everything is working... 
%   coarse nodes (cnds): doesn't seem like needs to be checked. 
%   coarse-2-fine (c2f) map: 
%       * It should map to only its subdomain 
%       * There should be no overlap of FEM nodes
%   Reg. Matrices (L):
%       * It subdomain should have no connections to any other
%       * The sum over rows and columns should equal zero
% 
%--------------------------------------------------------------------------
function check_softprior_cgrid_c2f_L(cnds,c2f,L,msh)
% c2f checks
% a. Get coarse grid indices
nsd = max(cnds(:,4));
figure
set(gcf,'position',[ 940         713        1328         420])
for n = 1:nsd
    subplot(1,nsd,n)
    is_sd{n} = find(cnds(:,4) == n);
    cvec           = zeros(size(cnds,1),1);
    cvec(is_sd{n}) = 1;
    fvec           = c2f*cvec;
    is_f{n}        = find( abs(fvec - 1) < 1e-3);
    plot3(msh.node(is_f{n},1),msh.node(is_f{n},2),msh.node(is_f{n},3),'.k','markersize',6)
end
subplot(1,nsd,1)
title('Associated FEM nodes for subdomain')
% b. check intersections
dins = zeros(nsd);
for n1 = 1:nsd
    for n2 = 1:nsd
        dins(n1,n2) = length(intersect(is_f{n1},is_f{n2}));
    end
end
disp('Check on intersection of FEM nodes between subdomains')
disp('Diagonals should equal number of FEM nodes per sudomain')
disp('Off-diagonals should be zero, i.e. no intersection')
dins
disp('The total number of FEM nodes should equal the sum of associated FEM nodes from each subdomain')
[size(msh.node,1) sum(diag(dins))]
%--------------------------------------------------------------------------
% L-regularization matrix checks
figure
vis_reg_connections(L,cnds,2,1)
title('Visualization of coarse node connections')
% b. checking sum over rows and columns
figure; hold on
plot(sum(L,1),'-k')
plot(sum(L,2),'--r')
title('Sum over rows and columns of L should be zero')





