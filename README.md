# S^2 Sampling Toolbox

Toolbox for generating uniform sampling patterns and decompositions of a unit sphere.

![S2 sampling demo](https://user-images.githubusercontent.com/13392426/57330923-6670c900-70e5-11e9-94e2-03756e51fd2e.jpg)

## Summary of Main Functions

>**`ParticleSampleSphere.m`**: generates an approximately uniform triangular tessellation of a unit sphere by using
 gradient descent to minimize a generalized electrostatic potential energy of a system of charged particles.
 In this implementation, initial configuration of particles is based on random sampling of a sphere, but 
 user-defined initializations are also permitted. Since the optimization algorithm implemented in this function 
 has O(N^2) complexity, it is not recommended that `ParticleSampleSphere.m` be used to optimize configurations 
 of more than 1E3 particles. Resolution of the meshes obtained with this function can be increased to an 
 arbitrary level with `SubdivideSphericalMesh.m`.

>**`SubdivideSphericalMesh.m`**: increases resolution of triangular or quadrilateral spherical meshes. Given a base
 mesh, its resolution is increased by a sequence of k subdivisions. Suppose that No is the original number of
 mesh vertices, then the total number of vertices after k subdivisions will be Nk=4^k*(No – 2)+2. This 
 relationship holds for both triangular and quadrilateral meshes.

>**`IcosahedronMesh.m`**: generates triangular surface mesh of an icosahedron. High-quality spherical meshes can be 
easily obtained by subdividing this base mesh with the `SubdivideSphericalMesh.m` function.

>**`QuadCubeMesh.m`**: generates quadrilateral mesh of a zero-centered unit cube. High-quality spherical meshes 
can be easily obtained by subdividing this base mesh with the `SubdivideSphericalMesh.m` function.

>**`SpiralSampleSphere.m`**: generates N uniformly distributed point samples on a unit sphere using a
 [spiral-based sampling method].

>**`RandSampleSphere.m`**: performs uniform random or stratified random sampling of a unit sphere with N point samples.

## Demo 1: Particle-based S^2 Sampling

	% Uniformly distribute 200 charged particles across the surface of a unit sphere.
	% This operation takes ~3 sec on my laptop (32GB RAM, 3.1 GHz CPU)
	[V,Tri,~,Ue]=ParticleSampleSphere('N',200);

	% Visualize optimization progress
	figure('color','w')
	plot(log10(1:numel(Ue)),Ue,'.-')
	set(get(gca,'Title'),'String','Optimization Progress','FontSize',40)
	set(gca,'FontSize',20)
	xlabel('log_{10}(Iteration #)','FontSize',30)
	ylabel('Electrostatic Potential','FontSize',30)

	% Visualize mesh based on computed configuration of particles
	figure('color','w')
	subplot(1,2,1)
	fv=struct('faces',Tri,'vertices',V);
	h=patch(fv);
	set(h,'EdgeColor','b','FaceColor','w')
	axis equal
	set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
	view(3)
	grid on
	set(get(gca,'Title'),'String','N=200 (base mesh)','FontSize',30)

	% Subdivide base mesh twice to obtain a spherical mesh of higher complexity
	fv_new=SubdivideSphericalMesh(fv,2);
	subplot(1,2,2)
	h=patch(fv_new);
	set(h,'EdgeColor','b','FaceColor','w')
	axis equal
	set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
	view(3)
	grid on
	set(get(gca,'Title'),'String','N=3170 (after 2 subdivisions)','FontSize',30)

## Demo 2: Icosahedron-based S^2 Decomposition

	% Get base icosahedron mesh
	TR=IcosahedronMesh; % also try this with `DodecahedronMesh.m`

	% Subdivide base mesh and visualize the results
	figure('color','w')
	ha=subplot(2,3,1);
	h=trimesh(TR); set(h,'EdgeColor','b','FaceColor','w')
	axis equal
	set(get(ha,'Title'),'String','base mesh (N=12)')
	for i=2:6
		ha=subplot(2,3,i);
		TR=SubdivideSphericalMesh(TR,1);
		h=trimesh(TR); set(h,'EdgeColor','b','FaceColor','w')
		axis equal
		set(get(ha,'Title'),'String',sprintf('N=%u',4^(i-1)*10+2))
		drawnow
	end

## Demo 3: Cuboid-based S^2 Decomposition

	% Get quad mesh of a unit cube
	fv=QuadCubeMesh; % also try this with `QuadRhombDodecMesh.m`

	% Subdivide the base mesh and visualize the results
	figure('color','w')
	ha=subplot(2,3,1);
	h=patch(fv); set(h,'EdgeColor','b','FaceColor','w')
	view(3)
	grid on
	axis equal
	set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
	set(get(ha,'Title'),'String','base mesh (N=8)')
	for i=2:6
		ha=subplot(2,3,i);
		fv=SubdivideSphericalMesh(fv,1);
		h=patch(fv); set(h,'EdgeColor','b','FaceColor','w')
		axis equal
		view(3)
		grid on
		set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
		set(get(ha,'Title'),'String',sprintf('N=%u',4^(i-1)*6+2))
		drawnow
	end
	
## License
[MIT] © 2019 Anton Semechko (a.semechko@gmail.com)


[spiral-based sampling method]: http://blog.wolfram.com/2011/07/28/how-i-made-wine-glasses-from-sunflowers/
[MIT]: https://github.com/AntonSemechko/S2-Sampling-Toolbox/blob/master/LICENSE.md