function TR=IcosahedronMesh
% Generate a triangular surface mesh of an icosahedron whose vertices lie
% on the surface of a zero-centered unit sphere.
%
%   - TR    : 'triangulation' object representing an icosahedron
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%

% Get the vertex coordinates
t=(1+sqrt(5))/2; % golden ratio
X=[0 1 t];
s=[1 1 1; 1 1 -1; 1 -1 -1; 1 -1 1];
X=repmat(X,[4 1]).*s;
X=[X;circshift(X,[0 -1]);circshift(X,[0 -2])];
x_L2=sqrt(sum(X.^2,2));
X=bsxfun(@rdivide,X,x_L2);

% Triangulate the points
Tri=convhull(X);
if ClosedMeshVolume({Tri X})<0, Tri=fliplr(Tri); end
TR=triangulation(Tri,X);
