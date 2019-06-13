function TR=IcosahedronMesh
% Generate a triangular surface mesh of an icosahedron whose vertices lie
% on the surface of a zero-centered unit sphere.
%
% OUTPUT:
%   - TR    : 'triangulation' object representing a zero-centered regular
%             icosahedron. This object has 12 evenly distributed vertices 
%             and 20 triangular faces.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Vertex coordinates
t=(1+sqrt(5))/2; % golden ratio
X=[0 1 t];
s=[1 1 1; 1 1 -1; 1 -1 -1; 1 -1 1];
X=repmat(X,[4 1]).*s;
X=[X;circshift(X,[0 -1]);circshift(X,[0 -2])];
X=ProjectOnSn(X);

% Triangulate the points
Tri=convhull(X);
if ClosedMeshVolume({Tri X})<0, Tri=fliplr(Tri); end
TR=triangulation(Tri,X);
