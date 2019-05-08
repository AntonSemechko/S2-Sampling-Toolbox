function [TR,W,G,SM,idx_unq]=TriQuad(TR,W,G)
% Subdivide a triangular surface mesh using generalized triangular 
% quadrisection. Triangular quadrisection is a linear subdivision procedure
% which inserts new vertices at the edge midpoints of the input mesh, 
% thereby producing four new faces for every face of the original mesh:
% 
%                     x3                        x3
%                    /  \      subdivision     /  \
%                   /    \        ====>       v3__v2
%                  /      \                  / \  / \
%                x1________x2              x1___v1___x2
%
%                      Original vertices : x1, x2, x3
%
%                      New vertices      : v1, v2, v3
%
%                      New faces         : [x1 v1 v3; x2 v2 v1; x3 v3 v2; v1 v2 v3] 
%
% In case of generalized triangular quadrisection, positions of the newly
% inserted vertices do not have to correspond to the edge midpoints, and 
% may be varied by assigning (positive) weights to the vertices of the 
% original mesh. For example, let xi and xj be two vertices connected by an
% edge, and suppose that Wi and Wj are the corresponding vertex weights. 
% Position of the new point on the edge (xi,xj) is defined as (Wi*xi+Wj*xj)/(Wi+Wj).
% Note that in order to avoid degeneracies and self-intersections, all 
% weights must be real numbers greater than zero.
%
% INPUT:
%   - TR   : surface mesh represented as an object of 'TriRep' class,
%            'triangulation' class, or a cell such that TR={Tri,X}, where
%            Tri is an M-by-3 array of faces and X is an N-by-3 array of 
%            vertex coordinates.
%   - W    : optional input argument specifying a N-by-1 array of STRICTLY
%            POSITIVE vertex weights used during interpolation of the new 
%            vertices, where N is the total number of the original mesh
%            vertices. 
%   - G    : optional input specifying scalar or vector field defined at 
%            the mesh vertices.
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as the input mesh.
%   - W   : interpolated vertex weights.
%   - G   : interpolated scalar or vector field defined at the vertices of 
%           the subdivided mesh.
%   - SM  : subdivision matrix.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Get the list of vertex co-ordinates and list of faces
[Tri,X,fmt]=GetMeshData(TR);

% Make sure that the mesh is composed entirely of triangles 
if size(Tri,2)~=3
    error('This function is meant for triangular surface meshes, not quad or tet meshes')
end

% Check vertex weights
if nargin<2 || isempty(W)
    W=[];
elseif ~ismatrix(W) || numel(W)~=size(X,1) || sum(W<=eps)>0
    error('W must be a N-by-1 array with positive entries, where N is the # of mesh vertices')
else
    W=W(:);
end

% Field?
if nargin==3 && ~isempty(G) && (size(G,1)~=size(X,1) || ~isnumeric(G) || ndims(G)>3)
    error('3-rd input argument must be a %u-by-d array where d>=1',size(X,1))
else
    G=[];
end


% Edges
E=[Tri(:,1) Tri(:,2); Tri(:,2) Tri(:,3); Tri(:,3) Tri(:,1)];
E=sort(E,2);
[E,idx_unq,idx]=unique(E,'rows','stable'); % setOrder='stable' ensures that identical results will be obtained for meshes with the same connectivity

% Compute new vertex positions
if ~isempty(W) % insert new vertices based on vertex weights
    
    w=bsxfun(@rdivide,[W(E(:,1)),W(E(:,2))],W(E(:,1))+W(E(:,2)));
    V=bsxfun(@times,X(E(:,1),:),w(:,1))+bsxfun(@times,X(E(:,2),:),w(:,2));
    
    if ~isempty(G) && nargout>2
        G=cat(1,G,bsxfun(@times,G(E(:,1),:,:),w(:,1))+bsxfun(@times,G(E(:,2),:,:),w(:,2)));
    end
    
    if nargout>1
        W=[W;W(E(:,1)).*w(:,1)+W(E(:,2)).*w(:,2)];
    end
        
else % insert new vertices at the edge mid-points
    
    V=(X(E(:,1),:)+X(E(:,2),:))/2;   
    if ~isempty(G) && nargout>2
        G=cat(1,G,(G(E(:,1),:,:)+G(E(:,2),:,:))/2);
    end
    
end

% Generate a subdivision matrix if one is required
Nx=size(X,1);   % # of vertices
Nt=size(Tri,1); % # of faces
if nargout>3
    
    dNx=size(E,1);
    if isempty(W), w=repmat([1 1]/2,[dNx 1]); end

    i=(1:dNx)'+Nx;
    i=cat(1,(1:Nx)',i,i);    
    j=cat(1,(1:Nx)',E(:));    
    w=cat(1,ones(Nx,1),w(:));
    
    SM=sparse(i,j,w,Nx+dNx,Nx,2*dNx+Nx);    
    
end

% Assign indices to new triangle vertices
V1= Nx + idx(1:Nt);
V2= Nx + idx((Nt+1):2*Nt);
V3= Nx + idx((2*Nt+1):3*Nt);

% Connectivities of the new faces
T1= [Tri(:,1) V1 V3];
T2= [Tri(:,2) V2 V1];
T3= [Tri(:,3) V3 V2];
T4= [V1       V2 V3];

T1=permute(T1,[3 1 2]);
T2=permute(T2,[3 1 2]);
T3=permute(T3,[3 1 2]);
T4=permute(T4,[3 1 2]);

Tri=cat(1,T1,T2,T3,T4);
Tri=reshape(Tri,[],3,1);

% New mesh
X=[X;V]; 
switch fmt
    case 1
        TR=triangulation(Tri,X);
    case 2
        TR=TriRep(Tri,X); %#ok<*DTRIREP>
    case 3
        TR={Tri X};
    case 4
        TR=struct('faces',Tri,'vertices',X);
end

