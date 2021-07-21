%--------------------------------------------------------------------------
%
% Merge coarse grids, maps, and regularization matrices
%
%--------------------------------------------------------------------------
function [cnds,c2f,L] = merge_grids_maps_Lmats(cgrids,c2fs,Ls)


%--------------------------------------------------------------------------
if length(cgrids) == 1
    error('There''s nothing to merge!')
elseif length(cgrids) == 2
    L    = [Ls{1} zeros( size(Ls{1},1), size(Ls{2},2)); ...
        zeros(size(Ls{2},1),size(Ls{1},2)) Ls{2}];        
    c2f  = [c2fs{1} c2fs{2}];
    % In addition to combined coarse grid node, add a flag for each
    % subdomain
    cnds = [cgrids{1}.node 1+0*cgrids{1}.node(:,1); ...
        cgrids{2}.node 2+0*cgrids{2}.node(:,1)];        
elseif length(cgrids) == 3
    L    = [Ls{1} zeros(size(Ls{1},1),size(Ls{2},2)) zeros(size(Ls{1},1),size(Ls{3},2)); ...
        zeros(size(Ls{2},1),size(Ls{1},2)) Ls{2} zeros(size(Ls{2},1),size(Ls{3},2)); ...
        zeros(size(Ls{3},1),size(Ls{1},2)) zeros(size(Ls{3},1),size(Ls{2},2)) Ls{3}];
    c2f  = [c2fs{1} c2fs{2} c2fs{3}];
    % In addition to combined coarse grid node, add a flag for each
    % subdomain
    cnds = [cgrids{1}.node 1+0*cgrids{1}.node(:,1); ...
        cgrids{2}.node 2+0*cgrids{2}.node(:,1); ...
        cgrids{3}.node 3+0*cgrids{3}.node(:,1)];
else
    error('not coded... but as long as the length isn''t to big, its easy to hardcode')
end
