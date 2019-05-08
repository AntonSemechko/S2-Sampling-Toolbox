function [fv,W,G]=QuadQuad(fv,W,G)
% Subdivide a quadrilateral surface mesh using generalized quadrilateral 
% quadrisection. This procedure is similar to generalized triangular 
% quadrisection described in the documentation of 'TriQuad' function, the 
% only difference being is that it uses quad surface meshes as input.
%
% INPUT:
%   - fv  :  input mesh specified as a face-vertex structure with fields 
%            'vertices' and 'faces' so that fv.vertices contains a N-by-3 
%            list of vertex co-ordinates and fv.faces contains a M-by-4 
%            list of faces. Alternatively, fv can be specified as a cell so 
%            that fv={F X}, where X is a N-by-3 list of vertex co-ordinates 
%            and F is a M-by-4 list of faces.
%   - W    : optional input argument. N-by-1 array of NON-ZERO, POSITIVE 
%            vertex weights used during interpolation of the new vertices,
%            where N is the total number of the original mesh vertices. 
%   - G    : scalar or vector field defined on the mesh vertices (optional).
%
% OUTPUT:
%   - fv   : subdivided mesh.
%   - W    : new set of vertex weights.   
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Get the list of vertex co-ordinates and list of faces
fmt=true;
if isstruct(fv) && sum(isfield(fv,{'vertices' 'faces'}))==2
    F=fv.faces;
    X=fv.vertices;
elseif iscell(fv) && numel(fv)==2
    F=fv{1};
    X=fv{2};      
    fmt=false;
else
    error('Unrecognized input format')
end

% Make sure that the mesh is composed entirely of quads 
if size(F,2)~=4 || ~isequal(round(F),F)
    error('Invalid entry for the list of faces')
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
E=[F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,4); F(:,4) F(:,1)];
E=sort(E,2);
[E,~,idx]=unique(E,'rows','stable'); % setOrder='stable' ensures that identical results will be obtained for meshes with the same connectivity

% Compute new vertex positions
if ~isempty(W)

    w=bsxfun(@rdivide,[W(E(:,1)),W(E(:,2))],W(E(:,1))+W(E(:,2)));
    w5=W(F); w5=bsxfun(@rdivide,w5,sum(w5,2));   
    
    V5=bsxfun(@times,w5(:,1),X(F(:,1),:)) + ...
       bsxfun(@times,w5(:,2),X(F(:,2),:)) + ...
       bsxfun(@times,w5(:,3),X(F(:,3),:)) + ...
       bsxfun(@times,w5(:,4),X(F(:,4),:));
    
    V=cat(1,bsxfun(@times,X(E(:,1),:),w(:,1))+bsxfun(@times,X(E(:,2),:),w(:,2)),V5);
   
    if ~isempty(G) && nargout>2
        G5=bsxfun(@times,w5(:,1),G(F(:,1),:,:)) + ...
           bsxfun(@times,w5(:,2),G(F(:,2),:,:)) + ...
           bsxfun(@times,w5(:,3),G(F(:,3),:,:)) + ...
           bsxfun(@times,w5(:,4),G(F(:,4),:,:));
        G=cat(1,G,bsxfun(@times,G(E(:,1),:,:),w(:,1))+bsxfun(@times,G(E(:,2),:,:),w(:,2)),G5);
    end    

else

    V5=(X(F(:,1),:)+X(F(:,2),:)+X(F(:,3),:)+X(F(:,4),:))/4;    
    V=cat(1,(X(E(:,1),:)+X(E(:,2),:))/2,V5);

    if ~isempty(G) && nargout>2
        G5=(G(F(:,1),:,:)+G(F(:,2),:,:)+G(F(:,3),:,:)+G(F(:,4),:,:))/4;
        G=cat(1,G,(G(E(:,1),:,:)+G(E(:,2),:,:))/2,G5);
    end
end

% Assign indices to the new vertices
Nx=size(X,1);   % # of vertices
Nt=size(F,1);   % # of faces
Ne=size(E,1);   % # of edges    

V1= Nx + idx(1:Nt);
V2= Nx + idx((Nt+1):2*Nt);
V3= Nx + idx((2*Nt+1):3*Nt);
V4= Nx + idx((3*Nt+1):4*Nt);
V5= Nx + Ne + (1:Nt)';

% Connectivities of the new faces
T1= [F(:,1) V1 V5 V4];
T2= [V1 F(:,2) V2 V5];
T3= [V5 V2 F(:,3) V3];
T4= [V4 V5 V3 F(:,4)];

T1=permute(T1,[3 1 2]);
T2=permute(T2,[3 1 2]);
T3=permute(T3,[3 1 2]);
T4=permute(T4,[3 1 2]);

F=cat(1,T1,T2,T3,T4);
F=reshape(F,[],4,1);

% New mesh
if fmt
    clear fv
    fv.faces=F;
    fv.vertices=[X;V];
else
    fv={F [X;V]};
end

