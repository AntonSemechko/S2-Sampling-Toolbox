function [V,Tri,Ue_i,Ue]=ParticleSampleSphere(varargin)
% Generate an approximately uniform triangular tessellation of a unit 
% sphere by minimizing generalized electrostatic potential energy 
% (i.e., Reisz s-energy) of a system of charged particles. Effectively, 
% this function produces a locally optimal solution to the problem that 
% involves finding a minimum Reisz s-energy configuration of N equally 
% charged  particles confined to the surface of a unit sphere; s=1 
% corresponds to the problem originally posed by J.J. Thomson. 
%
% SYNTAX:
% [V,Tri,Ue_i,Ue]=ParticleSampleSphere(option_i,value_i,option_j,value_j,...)
%
% OPTIONS:
%   - 'N'    : desired number of particles. Corresponding value of N must 
%              be a positive integer greater than 9. N=200 particles is the
%              default setting. A warning will be issued if N exceeds 1E3. 
%   - 'Vo'   : particle positions used to initialize the search.
%              Corresponding value of Vo must be a N-by-3 array, where N is
%              the number of particles. N=10 is the lowest permissible 
%              number of particles. Initializations consisting of more than 
%              1E3 particles are admissible, but may lead to unreasonably 
%              long optimization times.
%   - 's'    : Reisz s-energy parameter used to control the strength of
%              particle interactions. s must be a real number greater than 
%              zero. s=1 is the default setting.
%   - 'Etol' : absolute energy convergence tolerance. Etol=1E-5 is the 
%              default setting. Optimization will terminate when change in
%              potential energy between two consecutive iterations falls 
%              below Etol.
%   - 'Dtol' : maximum particle displacement tolerance (in degrees).
%              Dtol=1E-4 is the default setting. Optimization will
%              terminate when maximum displacement of any particle
%              between two consecutive iteration is less than Dtol.
%   - 'Nitr' : maximum number of iterations. Nitr must be a non-negative 
%              integer. Nitr=1E4 is the default setting. 
%
%              REMAINS TO BE IMPLEMENTED 
%   - 'CO'   : connectivity optimization setting. To obtain particle 
%              configurations with fewer dislocations (i.e., vertices      
%              possessing either less or more than 6 neighbours) set CO to
%              true. CO=false is the default setting.
%
% OUTPUT:
%   - V     : N-by-3 array of vertex (i.e., particle) coordinates.
%   - Tri   : M-by-3 list of face-vertex connectivities. 
%   - Ue_i  : N-by-1 array of particle energies.
%   - Ue    : k-by-1 array of potential energy values, where k-1 is the 
%             total number of iterations. Ue(1) and Ue(k) correspond to the
%             potential energies of the initial and final particle 
%             configurations, respectively.
%
% EXAMPLE: 
% -------------------------------------------------------------------------
% % Uniformly distribute 100 particles across the surface of a unit sphere
%
% [V,Tri,~,Ue]=ParticleSampleSphere('N',100);
% figure('color','w')
% subplot(1,2,1)
% plot(0:numel(Ue)-1,Ue,'.-')
% set(get(gca,'Title'),'String','Optimization Progress','FontSize',20)
% xlabel('Iteration #','FontSize',15)
% ylabel('Reisz s-enrgy','FontSize',15)
% subplot(1,2,2)
% h=patch('faces',Tri,'vertices',V); 
% set(h,'EdgeColor','b','FaceColor','w')
% axis equal vis3d
% view(3)
% grid on
% set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
% set(get(gca,'Title'),'String','Final Mesh','FontSize',20)
% -------------------------------------------------------------------------
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Check the inputs
[V,s,Etol,Dtol,Nitr,CO]=VerifyInputArgs(varargin); %#ok<*ASGLU>
if Nitr<0
    [Tri,Ue_i,Ue]=deal([]);
    return
end
N=size(V,1);    % number of particles
Dtol=Dtol/180*pi;
clear varargin
    
% Compute geodesic distances between particle pairs
DOT=V*V';       % dot product 
DOT(DOT<-1)=-1; DOT(DOT>1)=1;
GD=acos(DOT);   % geodesic distance

% Evaluate potential energy
GD(1:(N+1):end)=Inf; % set diagonal entries to Inf
Ue_ij=1./((GD.^s)+eps);
Ue_i=sum(Ue_ij,2);
Ue=sum(Ue_i);

% Approximate average distance between two neighbouring particles
d_ave=sqrt(8/sqrt(3)*pi/(size(V,1)-2));
if acos(1-d_ave^2/2)<1E-12
    fprintf(2,'Number of particles exceeds critical limit. Unable to continue due to limited precision of ''acos'' function.\n')
    [Tri,Ue_i,Ue]=deal([]);
    return
end
d_thr=max(1E-3*d_ave,1E-14);

% Iteratively optimize particle positions along negative gradient of
% potential energy using an adaptive Gauss-Seidel update scheme 
% -------------------------------------------------------------------------

fprintf('\nWait while particle positions are being optimized ...\n')
t0=clock;
i=0; 

a_min=1E-15;           % minimum step size
a_max=1;               % maximum step size
a=a_max*ones(N,1)/2;   % step sizes used during position updates

idx_jo=true(N,1);
[dE,dV]=deal(Inf);
Vrec=repmat({V},[1 10]);
while i<Nitr && dE>Etol && dV>Dtol
    
    i=i+1;
    
    % Sort particles according to their energy contribution
    [~,idx_sort]=sort(abs(Ue_i-mean(Ue_i)),'descend');
     
    % Update positions of individual particles
    dV_max=0;
    for k=1:N
        
        j=idx_sort(k);
        
        idx_j=idx_jo; 
        idx_j(j)=false; % particle indices, except the current one 
        
        % Potential energy gradient of the j-th particle
        DOTj=DOT(idx_j,j);
        GDj=GD(idx_j,j);
        
        dVj=bsxfun(@times,s./(sqrt(1-DOTj.^2)+eps),V(idx_j,:));                
        dVj=bsxfun(@rdivide,dVj,(GDj.^(s+1))+eps);
        
         if min(GDj)<d_thr
             
            % remove influence of particles closer than d_thr
            idx_fail=GDj<d_thr;
            dVj(idx_fail,:)=[];
            
            if isempty(dVj) 
                
                % perturb positions of all particles
                V=V+randn(size(V))/1E3;
                V=ProjectOnSn(V);
                
                DOT=V*V';       
                DOT(DOT<-1)=-1; DOT(DOT>1)=1;
                GD=acos(DOT);   
                GD(1:(N+1):end)=Inf; 
                Ue_ij=1./((GD.^s)+eps);
                Ue_i=sum(Ue_ij,2);

                dV_max=pi;
                break
            end 
            
         end        
        dVj=sum(dVj,1);
        
        % Only retain tangential component of the gradient
        dVj_n=(dVj*V(j,:)')*V(j,:);
        dVj_t=dVj-dVj_n;        

        % Update position of the j-th particle
        m=0; 
        Uj_old=sum(Ue_ij(j,:));
        while m<50
        
            m=m+1;
            
            % Update position of the j-th particle
            Vj_new=V(j,:)-a(j)*dVj_t;
            Vj_new=Vj_new/norm(Vj_new);
                     
            % Recompute dot products and geodesic distances
            DOTj=sum(bsxfun(@times,V,Vj_new),2);
            DOTj(DOTj<-1)=-1; DOTj(DOTj>1)=1;
            GDj=acos(DOTj);
            GDj(j)=Inf;
            
            Ue_ij_j=1./((GDj.^s)+eps);
            
            % Check if the system potential decreased
            if sum(Ue_ij_j)<Uj_old
                
                dV_max=max(dV_max,norm(Vj_new-Vrec{1}(j,:)));
                V(j,:)=Vj_new;
                
                DOT(j,:)=DOTj';
                DOT(:,j)=DOTj;
                
                GD(j,:)=GDj';
                GD(:,j)=GDj;
                
                Ue_ij(j,:)=Ue_ij_j';
                Ue_ij(:,j)=Ue_ij_j;
                
                if m==1, a(j)=min(a(j)*1.05,a_max); end
                
                break
                
            else
                
                if a(j)>a_min
                    a(j)=max(a(j)/1.1,a_min);
                else
                    break
                end
                
            end
            
        end
        
    end
    
    % Evaluate net potential energy of the system 
    Ue_i=sum(Ue_ij,2);
    Ue(i+1)=sum(Ue_i);
    
    % Average change in potential energy
    if i>=10, dE=(Ue(end-10)-Ue(end))/10; end
    
    % Maximum displacement (in radians)
    dV_max=0.1*acos(1-dV_max^2/2); 
    dV=dV_max;
    Vrec{10}=V;
    Vrec=circshift(Vrec,[0 -1]);
    
    % Reset the step sizes; to avoid premature convergence
    if mod(i,40)==0
        a_max=min(a_max,5*max(a));
        a(:)=a_max;
    end

    % Progress update
    if (mod(i,400)==0 && i>100) || i==1
        fprintf('\n%-15s %-15s     %-15s     %-15s     %-15s\n','Iteration #','log(Energy/N)','log(dE/Etol)','log(dV/Dtol)','Time (sec)')
    end
    
    if mod(i,10)==0 || i==1
        fprintf('%-15u %-15.11f     %-15.2f     %-15.2f     %-15.1f\n',i,log10(Ue(end)/size(V,1)),log10(dE/Etol),log10(dV_max/Dtol),etime(clock,t0))
    end

end
clear DOT GD Ue_ij


if mod(i,10)~=0 && i>1
    fprintf('%-15u %-15.11f     %-15.2f     %-15.2f     %-15.1f\n',i,log10(Ue(end)/size(V,1)),log10(dE/Etol),log10(dV_max/Dtol),etime(clock,t0))
end

if Nitr>0
    fprintf('Optimization terminated after %u iterations. Elapsed time: %5.1f sec\n',i,etime(clock,t0))
    fprintf('Convergence thresholds: Etol=%.6e, Dtol=%.6e degrees\n',Etol,Dtol/pi*180)
end

% Triangulate particle positions
if nargout>1    
    Tri=convhull(V);
    if ClosedMeshVolume({Tri V})<0, Tri=fliplr(Tri); end
end


%==========================================================================
function [V,s,Etol,Dtol,Nitr,CO]=VerifyInputArgs(VarsIn)
% Make sure user-defined input arguments have valid format

% Default settings
V=[];
s=1; 
Etol=1E-5; 
Dtol=1E-4; 
Nitr=1E4;
CO=false;
if isempty(VarsIn)
    V=RandSampleSphere;
    return; 
end

% Check that there is an even number of inputs
Narg=numel(VarsIn);
if mod(Narg,2)~=0
    error('Incorrect number of input arguments')
end

% Check user-defined parameters
FNo={'N','Vo','s','Etol','Dtol','Nitr','Alg'};
flag=false(1,7); exit_flag=false;
for i=1:Narg/2
    
    % Make sure the input is a string
    str=VarsIn{2*(i-1)+1};
    if ~ischar(str)
        error('Input argument #%u is not a valid paramer name',2*(i-1)+1)
    end
    
    % Get parameter "value"
    Val=VarsIn{2*i};
    
    % Match the string against the list of available options  
    chk=strcmpi(str,FNo);
    id=find(chk,1);
    if isempty(id), id=0; end
    
    switch id
        case 1 % number of particles
            
            % Check if 'initialization' option has also been specified
            if flag(2) 
                error('Ambiguous combination of options. Specify option ''%s'' or option ''%s'', but not both.',FNo{2},FNo{1})
            end
            
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val)|| Val<10 || ~isfinite(Val) || ~isreal(Val)   
                error('Incorrect entry for the ''%s'' option. N must be a positive integer greater than 9.',FNo{1})
            end
            V=RandSampleSphere(round(Val)); %#ok<*NASGU>
            
            % Check if there are more than 1E3 particles
            if size(V,1)>1E3
                
                % Construct a 'yes'/'no' questdlg
                choice = questdlg('Default particle limit exceeded. Would you like to continue?', ...
                    'Particle Limit Exceeded','   YES   ','   NO   ','   NO   ');
                
                % Handle response
                if strcmpi(choice,'   NO   '), exit_flag=true; end
                
            end
            
        case 2 % initialization
        
            % Check if 'number' option has also been specified
            if flag(1) 
                error('Ambiguous combination of options. Specify option ''%s'' or option ''%s'', but not both.',FNo{1},FNo{2})
            end
            
            % Check the format
            if ~ismatrix(Val) || ~isnumeric(Val) || size(Val,2)~=3 || size(Val,1)<10 || sum(~isfinite(Val(:)))>0 || ~isreal(Val)
                error('Incorrect entry for the ''%s'' option. Vo must be a N-by-3 array, where N>9 is the number of particles.',FNo{2})
            end
                        
            % Make sure the particles are on the unit sphere
            V=ProjectOnSn(Val);
            clear Val V_L2 
            
            % Check if there are more than 1E3 particles
            if size(V,1)>1E3
                
                % Construct a 'yes'/'no' questdlg
                choice = questdlg('Default particle limit exceeded. Would you like to continue?', ...
                    'Particle Limit Exceeded','   YES   ','   NO   ','   NO   ');
                
                % Handle response
                if strcmpi(choice,'   NO   '), exit_flag=true; end
                
            end
            
        case 3 % s parameter
            
            % Check the format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<1E-6 || ~isreal(Val)
                error('Incorrect entry for the ''%s'' option. s must be a positive real number.',FNo{3})
            end
            s=Val;
            
        case 4 % energy tolerance
        
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<eps || ~isreal(Val) 
                error('Incorrect entry for the ''%s'' option. Etol must be a positive real number.',FNo{4})
            end
            Etol=Val;
            
        case 5 % maximum dispalcement tolerance
            
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<eps || ~isreal(Val)
                error('Incorrect entry for the ''%s'' option. Dtol must be a positive real number.',FNo{5})
            end
            Dtol=Val; 
            
        case 6 % maximum number of iterations 
            
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<0 || ~isreal(Val)
                error('Incorrect entry for the ''%s'' option. Nitr must be a non-negative integer.',FNo{6})
            end
            Nitr=Val;
            
        case 7 % connectivity optimization
            
            % Check format
            if numel(Val)~=1 || ~islogical(Val)
                error('Incorrect entry for the ''%s'' option. OC must be either true or false.',FNo{7})
            end
            CO=Val;
            
        otherwise
            error('''%s'' is not a recognized parameter name',str)
    end
    flag(id)=true;
    
end
if isempty(V), V=RandSampleSphere; end

if exit_flag, Nitr=-1; end %#ok<*UNRCH>

