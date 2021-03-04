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
% INPUT PARAMETERS/OPTIONS:
%   - 'N'    : positive integer specifying either the desired number of 
%              particles or antipodal particle pairs ( see 'asym' option 
%              below). Default settings for 'N' is 200 when 'asym'=false 
%              and 100 when 'asym'=true. If 'N' exceeds 1E3, user will be 
%              prompted to manually verify if he wishes to continue. Set 
%              'qdlg' to false to disable manual verification when N>1E3. 
%              The lowest permissible number of particles is 14.  
%   - 'Vo'   : array of particle positions used to initialize the search.
%              This input is optional and should be used in place of 'N' 
%              when suboptimal initial configuration of particles is 
%              available. Corresponding value of 'Vo' must be a N-by-3 
%              array of particle coordinates, where N the number of 
%              particles (or antipodal particle pairs if 'asym' is true). 
%              Note that initializations consisting of more than 1E3 
%              particles (or particle pairs) are admissible, but may lead 
%              to unreasonably long optimization times.
%   - 's'    : Reisz s-energy parameter used to control the strength of
%              particle interactions; higher values of 's' lead to stronger
%              short-range interactions . 's' must be a real number greater
%              than zero. 's'=1 is the default setting.
%   - 'asym' : compute antipodally symmetric particle configurations. Set 
%              'asym' to true to obtain a uniformly distributed set of 2N 
%              particles comprised of N antipodal particle pairs. Recall, 
%              an antipodal partner of particle P is -P (i.e., P reflected 
%              through the origin). 'asym'=false is the default setting. 
%   - 'Etol' : absolute energy convergence tolerance. Optimization will 
%              terminate when change in potential energy between two 
%              consecutive iterations falls below Etol. 'Etol'=1E-5 is the 
%              default setting.
%   - 'Dtol' : maximum particle displacement tolerance (in degrees).
%              Optimization will terminate when maximum displacement of any
%              particle between two consecutive iteration is less than 
%              Dtol. 'Dtol'=1E-4 is the default setting.
%   - 'Nitr' : maximum number of iterations. Nitr must be a non-negative 
%              integer. 'Nitr'=1E4 is the default setting. 
%   - 'upd'  : progress update.  Set 'upd' to false to disable progress
%              updates. 'upd'=true is the default setting.
%   - 'qdlg' : default maximum particle limit verification. Set 'qdlg' to 
%              false to disable the question dialog pop-up prompting the 
%              user to indicate if they wish to continue when N>1E3.
%              'qdlg'=true is the default setting.
%
%              REMAINS TO BE IMPLEMENTED 
%   - 'CO'   : connectivity optimization. To obtain particle configurations
%              with fewer dislocations (i.e., vertices possessing either 
%              less or more than 6 neighbours) set 'CO' to true. 'CO'=false
%              is the default setting.
%
% OUTPUT:
%   - V     : N-by-3 array of vertex (i.e., particle) coordinates. When
%             'asym'=true, -V(i,:) is the antipodal partner of V(i,:). 
%   - Tri   : M-by-3 list of face-vertex connectivities. When 'asym'=false,
%             Tri is triangulation of V. When 'asym'=true, Tri is
%             triangulation of 2N particles [V;-V].
%   - Ue_i  : N-by-1 array of particle (or particle pair) energies.
%   - Ue    : k-by-1 array of potential energy values, where k-1 is the 
%             total number of iterations. Ue(1) and Ue(k) correspond to the
%             potential energies of the initial and final particle 
%             configurations, respectively.
%
%
% EXAMPLE 1: Uniformly distribute 200 particles across the surface of a 
%            unit sphere
% -------------------------------------------------------------------------
% % Sample
% [V,Tri,~,Ue]=ParticleSampleSphere('N',200);
%
% % Visualize optimization progress
% figure('color','w')
% subplot(1,2,1)
% plot(log10(1:numel(Ue)),Ue,'.-')
% set(get(gca,'Title'),'String','Optimization Progress','FontSize',40)
% set(gca,'FontSize',20,'XColor','k','YColor','k')
% xlabel('log_{10}(Iteration #)','FontSize',30,'Color','k')
% ylabel('Reisz s-Energy','FontSize',30,'Color','k')
%
% % Visualize mesh
% subplot(1,2,2)
% h=patch('faces',Tri,'vertices',V);
% set(h,'EdgeColor','b','FaceColor','w')
% axis equal
% hold on
% plot3(V(:,1),V(:,2),V(:,3),'.k','MarkerSize',15)
% set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
% view(3)
% grid off
% set(get(gca,'Title'),'String','N=200 (base mesh)','FontSize',30)
% -------------------------------------------------------------------------
%
%
% EXAMPLE 2: Uniformly distribute 100 antipodally symmetric particle pairs
%            across the surface of a unit sphere
% -------------------------------------------------------------------------
% %Sample
%[V,Tri,~,Ue]=ParticleSampleSphere('N',100,'asym',true);

% %Visualize optimization progress
%figure('color','w')
%subplot(1,2,1)
%plot(log10(1:numel(Ue)),Ue,'.-')
%set(get(gca,'Title'),'String','Optimization Progress','FontSize',40)
%set(gca,'FontSize',20,'XColor','k','YColor','k')
%xlabel('log_{10}(Iteration #)','FontSize',30,'Color','k')
%ylabel('Reisz s-Energy','FontSize',30,'Color','k')

% %Visualize mesh. Note that vertices of the mesh are [V;-V] and not [V]
% %used in the previous example. This is because -V are antipodal partners
% %of V. However, just like in the previous example, computed mesh is also
% composed of 200 vertices. 
%subplot(1,2,2)
%fv=struct('faces',Tri,'vertices',[V;-V]); 
%h=patch(fv);
%set(h,'EdgeColor','b','FaceColor','w')
%axis equal
%set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
%view(3)
%grid off
%hold on
%plot3(V(:,1),V(:,2),V(:,3),'.k','MarkerSize',15)
%plot3(-V(:,1),-V(:,2),-V(:,3),'.r','MarkerSize',15)
%set(get(gca,'Title'),'String','Final Mesh','FontSize',30)
% -------------------------------------------------------------------------
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Check the inputs
prms=VerifyInputArgs(varargin); %#ok<*ASGLU>
[V,s,upd,CO]=deal(prms.V,prms.s,prms.upd,prms.CO);

if prms.Nitr<0
    [Tri,Ue_i,Ue]=deal([]);
    return
end
N=size(V,1);    % number of particles
prms.Dtol=prms.Dtol/180*pi;
clear varargin

if prms.Nitr>=0 || nargout>2
    
    % Compute geodesic distances between particle pairs  
    DOT=max(min(V*V',1),-1); % dot product
    GD=acos(DOT);            % geodesic distance
    
    % Evaluate potential energy
    GD(1:(N+1):end)=Inf; % set diagonal entries to Inf
    Ue_ij=1./(eps + GD.^s);
    if prms.asym
        Ue_ij=Ue_ij + 1./(eps + (pi-GD).^s);        
    end
    Ue_i=sum(Ue_ij,2);
    if prms.asym
        Ue_i=Ue_i+(1/pi)^s;
    end
    Ue=sum(Ue_i);
    
    % Approximate average distance between two neighbouring particles
    if prms.asym
        d_ave=sqrt(8/sqrt(3)*pi/(2*N-2));
    else
        d_ave=sqrt(8/sqrt(3)*pi/(N-2));
    end    
    if acos(1-d_ave^2/2)<1E-12
        fprintf(2,'Number of particles exceeds critical limit. Unable to continue due to limited numerical precision of ''acos'' function.\n')
        [Tri,Ue_i,Ue]=deal([]);
        return
    end
    d_thr=max(1E-3*d_ave,1E-14);
    
end

% Iteratively optimize particle positions along negative gradient of
% potential energy using an adaptive Gauss-Seidel update scheme 
% -------------------------------------------------------------------------
if upd && prms.Nitr>0
    fprintf('\nWait while particle positions are being optimized ...\n')
    t0=clock;
end

a_min=1E-15;           % minimum step size
a_max=1;               % maximum step size
a=a_max*ones(N,1)/2;   % step sizes used during position updates

idx_jo=true(N,1);
[dE,dV]=deal(Inf);
Vrec=repmat({V},[1 10]);

i=0;
while i<prms.Nitr && dE>prms.Etol && dV>prms.Dtol
    
    i=i+1;
    
    % Sort particles according to their energy contribution
    [~,idx_sort]=sort(abs(Ue_i-mean(Ue_i)),'descend');
     
    % Update positions of individual particles (or particle pairs)
    dV_max=0;
    for k=1:N
        
        j=idx_sort(k);
        
        idx_j=idx_jo; 
        idx_j(j)=false; % particle indices, except the current one 
        
        % Potential energy gradient of the j-th particle
        DOTj=DOT(idx_j,j);
        GDj=GD(idx_j,j);
        
        dVj=bsxfun(@times,s./(eps + sqrt(1-DOTj.^2)),V(idx_j,:));
        if prms.asym
            dVj=bsxfun(@rdivide,dVj,eps + GDj.^(s+1)) - bsxfun(@rdivide,dVj,eps + (pi-GDj).^(s+1));
        else
            dVj=bsxfun(@rdivide,dVj,eps + GDj.^(s+1));
        end
        
        if min(GDj)<d_thr || (prms.asym && min(pi-GDj)<d_thr)
            
            % Remove influence of particles closer than d_thr
            idx_fail=GDj<d_thr;
            if prms.asym
                idx_fail=idx_fail | (pi-GDj)<d_thr;
            end
            dVj(idx_fail,:)=[];
            
            if isempty(dVj)
                
                % Perturb positions of all particles
                V=V+randn(size(V))/1E3;
                V=ProjectOnSn(V);
                
                DOT=max(min(V*V',1),-1);
                GD=acos(DOT);
                GD(1:(N+1):end)=Inf;
                Ue_ij=1./(eps + GD.^s);
                if prms.asym
                    Ue_ij=Ue_ij + 1./(eps + (pi-GD).^s);
                end
                Ue_i=sum(Ue_ij,2);
                if prms.asym
                    Ue_i=Ue_i+(1/pi)^s; %#ok<*NASGU>
                end
                
                dV_max=2;
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
            DOTj=max(min(V*Vj_new(:),1),-1);
            GDj=acos(DOTj);
            GDj(j)=Inf;

            Ue_ij_j=1./(eps + GDj.^s);
            if prms.asym
                Ue_ij_j=Ue_ij_j + 1./(eps + (pi-GDj).^s);
            end
            
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
    if prms.asym, Ue_i=Ue_i + (1/pi)^s; end
    Ue(i+1)=sum(Ue_i);
    
    % Maximum displacement (in radians)
    dV_max=acos(1-dV_max^2/2); 
    Vrec{10}=V;
    Vrec=circshift(Vrec,[0 -1]);
    
    % Average change in potential energy
    if i>=10
        [dE,dE2]=deal((Ue(end-10)-Ue(end))/10);
        dV=dV_max;
    elseif i==1
        dE2=(Ue(1)-Ue(end))/i;
    end
    
    % Progress update
    if upd && ((mod(i,400)==0 && i>100) || i==1)
        fprintf('\n%-15s %-15s     %-15s     %-15s     %-15s\n','Iteration #','log(Energy/N)','log(dE/Etol)','log(dV/Dtol)','Time (sec)')
    end
    
    if upd && (mod(i,10)==0 || i==1)
        fprintf('%-15u %-15.11f     %-15.2f     %-15.2f     %-15.1f\n',i,log10(Ue(end)/size(V,1)),log10(max(dE2,eps)/prms.Etol),log10(max(dV_max,eps)/prms.Dtol),etime(clock,t0))
    end
    
    % Reset step sizes; to avoid premature convergence
    if mod(i,40)==0
        a_max=min(a_max,5*max(a));
        a(:)=a_max;
    end
    
end
clear DOT GD Ue_ij

if upd && (mod(i,10)~=0 && i>1)
    fprintf('%-15u %-15.11f     %-15.2f     %-15.2f     %-15.1f\n',i,log10(Ue(end)/size(V,1)),log10(dE/prms.Etol),log10(dV_max/prms.Dtol),etime(clock,t0))
end

if upd && prms.Nitr>0
    fprintf('Optimization terminated after %u iterations. Elapsed time: %5.1f sec\n',i,etime(clock,t0))
    fprintf('Convergence tolerances: Etol=%.6e, Dtol=%.6e degrees\n',prms.Etol,prms.Dtol/pi*180)
end

if ~prms.user_init 
    if prms.asym
        idx=V(:,3)<0;
        V(idx,:)=-V(idx,:);
    end
    [~,id_srt]=sort(V(:,3),'descend');
    V=V(id_srt,:);
    Ue_i=Ue_i(id_srt);
end


% Triangulate particle positions
if nargout>1
    if prms.asym
        Tri=convhull([V;-V]);
        if ClosedMeshVolume({Tri [V;-V]})<0, Tri=fliplr(Tri); end        
    else
        Tri=convhull(V);
        if ClosedMeshVolume({Tri V})<0, Tri=fliplr(Tri); end
    end
end


%==========================================================================
function prms=VerifyInputArgs(VarsIn)
% Make sure user-defined input arguments have valid format

% Default settings
prms.V=[];
prms.s=1; 
prms.Etol=1E-5; 
prms.Dtol=1E-4; 
prms.Nitr=1E4;
prms.asym=false;
prms.upd=true;
prms.CO=false;
prms.qdlg=true;
prms.user_init=false;
if isempty(VarsIn)
    prms.V=RandSampleSphere;
    return; 
end

% Check that there is an even number of inputs
Narg=numel(VarsIn);
if mod(Narg,2)~=0
    error('This function only accepts name-value parameter pairs as inputs')
end

% Check user-defined parameters
FNo={'N','Vo','s','asym','Etol','Dtol','Nitr','upd','CO','qdlg'};
flag=false(1,numel(FNo)); exit_flag=false;
N=[];
for i=1:Narg/2
    
    % Make sure the input is a string
    str=VarsIn{2*(i-1)+1};
    if ~ischar(str) || numel(str)>4
        error('Input argument #%u is not a valid paramer name',2*(i-1)+1)
    end
    
    % Get parameter "value"
    Val=VarsIn{2*i};
    
    % Match the string against the list of available options  
    chk=strcmpi(str,FNo);
    id=find(chk,1);
    if isempty(id), id=0; end
    
    switch id
        case 1 % number of particles (or particle pairs)
            
            % Check if 'initialization' option has also been specified
            if flag(2) 
                error('Ambiguous combination of input parameters. Specify ''%s'' or ''%s'', but not both.',FNo{2},FNo{1})
            end
            N=Val;                        
            
        case 2 % initialization
        
            % Check if 'number' option has also been specified
            if flag(1) 
                error('Ambiguous combination of input parameters. Specify ''%s'' or ''%s'', but not both.',FNo{1},FNo{2})
            end
            
            % Check the format
            if ~ismatrix(Val) || ~isnumeric(Val) || size(Val,2)~=3 || size(Val,1)<14 || any(~isfinite(Val(:)))
                error('Incorrect entry for ''%s''. ''%s'' must be set to a N-by-3 array, where N is the number of particles (or particle pairs).',FNo{2},FNo{2})
            end
                        
            % Make sure particles are on the unit sphere
            prms.V=ProjectOnSn(Val);
            prms.user_init=true;
            N=size(prms.V,1);
            
        case 3 % s parameter
            
            % Check the format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<1E-6 
                error('Incorrect entry for the ''%s'' parameter. ''%s'' must be set to a positive real number.',FNo{3},FNo{3})
            end
            prms.s=Val;
            
        case 4 % antipodal particle pairs
            
            % Check format
            if numel(Val)~=1 || ~islogical(Val)
                error('Incorrect entry for ''%s''. ''%s'' must be set to true or false.',FNo{4},FNo{4})
            end
            prms.asym=Val;
            
        case 5 % energy tolerance
        
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<eps  
                error('Incorrect entry for the ''%s'' parameter. ''%s'' must be set to a positive real number.',FNo{5},FNo{5})
            end
            prms.Etol=Val;
            
        case 6 % maximum displacement tolerance
            
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<eps 
                error('Incorrect entry for the ''%s'' parameter. ''%s'' must be set to a positive real number.',FNo{6},FNo{6})
            end
            prms.Dtol=Val; 
            
        case 7 % maximum number of iterations 
            
            % Check format
            if numel(Val)~=1 || ~isnumeric(Val) || ~isfinite(Val) || Val<0 
                error('Incorrect entry for the ''%s'' parameter. ''%s'' must be set to a non-negative integer.',FNo{7},FNo{7})
            end
            prms.Nitr=Val;
            
        case 8 % progress update
            
            % Check format
            if numel(Val)~=1 || ~islogical(Val)
                error('Incorrect entry for ''%s''. ''%s'' must be set to true or false.',FNo{8},FNo{8})
            end
            prms.upd=Val;
            
        case 9 % connectivity optimization
            
            % Check format
            if numel(Val)~=1 || ~islogical(Val)
                error('Incorrect entry for ''%s''. ''%s'' must be set to true or false.',FNo{9},FNo{9})
            end
            prms.CO=Val;
            
        case 10 % disable question dialog pop-up when N>1E3
            
            % Check format
            if numel(Val)~=1 || ~islogical(Val)
                error('Incorrect entry for ''%s''. ''%s'' must be set to true or false.',FNo{10},FNo{10})
            end
            prms.qdlg=Val;
            
        otherwise
            error('''%s'' is not a recognized parameter name',str)
    end
    flag(id)=true;
    
end


if ~isempty(N)
    chk=numel(N)~=1 | ~isnumeric(N) | any(~isfinite(N(:))) | ~isequal(round(N),N);
    if  chk || (prms.asym && N<7) || (~prms.asym && N<14)
        error('Incorrect entry for ''N''. Total number of particles must be greater than 13.')
    end
else
    N=200;
end

if isempty(prms.V)
    prms.V=RandSampleSphere(N);
end

% Check if there are more than 1E3 particles
if N>1E3 && prms.qdlg
    
    % Construct a 'yes'/'no' questdlg
    choice = questdlg('Default particle limit exceeded. Would you like to continue?', ...
        'Particle Limit Exceeded','   YES   ','   NO   ','   NO   ');
    
    % Handle response
    if strcmpi(choice,'   NO   '), exit_flag=true; end
    
end

if exit_flag, prms.Nitr=-1; end %#ok<*UNRCH>

