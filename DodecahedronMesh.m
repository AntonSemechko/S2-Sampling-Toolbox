function TR=DodecahedronMesh
% Generate a triangular surface mesh of a Pentakis dodecahedron. This 
% object is a dual representation of an icosahedron; i.e., its vertices
% correspond to the centroids of icosahedron's faces.
%
%   - TR    : 'triangulation' object representing a Pentakis dodecahedron.           
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


tr=IcosahedronMesh;
X=tr.Points;
F=tr.ConnectivityList;

C=(X(F(:,1),:)+X(F(:,2),:)+X(F(:,3),:))/3;
C=ProjectOnSn(C);

X=[X;C];

Tri=convhull(X);
if ClosedMeshVolume({Tri X})<0, Tri=fliplr(Tri); end
TR=triangulation(Tri,X);

