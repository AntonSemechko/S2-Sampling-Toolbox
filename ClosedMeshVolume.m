function Vol=ClosedMeshVolume(TR)
% Compute volume of a region enclosed by a triangular surface mesh.
%
% INPUT:
%   - TR    : input surface mesh represented as an object of 'TriRep' 
%             class, 'triangulation' class, or a cell such that TR={Tri,V},
%             where Tri is an M-by-3 array of faces and V is an N-by-3 
%             array of vertex coordinates. 
%
% OUTPUT:
%   - Vol   : real number specifying volume enclosed by TR. Vol<0 means 
%             that the order of mesh vertices in the face connectivity list
%             is in the clockwise direction, as viewed from the outside the
%             mesh.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Face and vertex lists
[Tri,V]=GetMeshData(TR);
if size(Tri,2)~=3
    error('This function is meant ONLY for triangular surface meshes.')
end
    
% Face vertices
V1=V(Tri(:,1),:);
V2=V(Tri(:,2),:);
V3=V(Tri(:,3),:);

% Face centroids
C=(V1+V2+V3)/3;

% Face normals
FN=cross(V2-V1,V3-V1,2);

% Volume
Vol=sum(dot(C,FN,2))/6;

