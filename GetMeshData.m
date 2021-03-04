function [Tri,X,fmt]=GetMeshData(TR)
% Get face-vertex connectivity list and list of vertex coordinates of a
% mesh (TR) represented in one of the following formats:
%
%  (1) 'triangulation' object, 
%  (2) 'TriRep' object, 
%  (3) 1-by-2 cell such that TR={Tri,V}, where Tri is a M-by-3 or M-by-4
%      array of elements and V is a N-by-3 array of vertex coordinates, or
%  (4) structure with fields 'faces' and 'vertices', such that 
%      Tri=TR.faces, and V=TR.vertices, where Tri and V have the same
%      meaning as in (3).
%
% OUTPUT:
%   - Tri   : M-by-3 array of face-vertex connectivities when TR is 
%             triangular mesh, and M-by-4 array when TR is a quad or tet 
%             mesh.
%   - X     : N-by-3 array of vertex coordinates
%   - fmt   : mesh format (see above)
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Get face and vertex lists
if isa(TR,'triangulation')
    Tri=TR.ConnectivityList;
    X=TR.Points;
    fmt=1;
elseif isa(TR,'TriRep')
    Tri=TR.Triangulation;
    X=TR.X;
    fmt=2;
elseif iscell(TR) && numel(TR)==2
    Tri=TR{1};
    X=TR{2};
    fmt=3;
elseif isstruct(TR) && isfield(TR,'faces') && isfield(TR,'vertices')
    Tri=TR.faces;
    X=TR.vertices;
    fmt=4;
else
    error('Unrecognized mesh format')
end
if fmt<3, return; end

% Check that face and vertex lists have correct dimensions
c=size(Tri,2);
if ~isnumeric(Tri) || ~ismatrix(Tri) || c<3 || c>4 || ~isequal(Tri,round(Tri)) || any(Tri(:)<0)
    error('Vertex connectivity list must be specified as a M-by-(3 or 4) array of POSITIVE integers')
end

d=size(X,2);
if ~isnumeric(X) || ~ismatrix(X) || any(~isfinite(X(:))) || d<2 || d>3 
    error('List of vertex coordinates must be specified as a N-by-2 or N-by-3 array of real numbers')
end

