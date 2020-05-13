% The initialization and number of cores needed in the execution of the code
parpool(7);

% The different stresses in the phase diagram for DST suspensions 
stress = [0.1 0.5 1.0 2.0 5.0 10.0 20.0 100.0];

% The packing fractions in the phase diagram for DST suspensions
phi = [0.76 0.77 0.78 0.785 0.79 0.8 0.82];

% The number of particles that is in the simulation box
NN = [128 256 512 768 1024 1280 1408 1536];

% Use of different random numbers for initializing the seed for different parallel cores. 
rng('shuffle','twister');
rnd_numbs = (10^9).*rand(8,1)

% The box parameters for the initial configuration
mu = [1/3 1/3 1/3 1/3 1/3 1/3 1/3 1/3];
% THe increment in periodic vectors is decided by this variable
dstep = [1 10 20 100 200 500 1000];
% The temperature is decided by this variable
b1 = [1 1 1 1 1 1 1 1];

% for loop for going through different cores 
parfor pts1 = 1:1:7

    pts = 4;   
    % This array consists of pre-decided number of MC steps after which the 
    % relevant variables are saved into file.
    set = [1000:1000:(10^6)];
    
    % Useful dimensions required for the program
    lr = 50.0;
    lx = 50.0;
    ly = 50.0;
    lx1 = 10.0;
    ly1 = 10.0;
    
    dhx = 0.1;
    dhy = 0.1;
    
    dh = 0.001;
    
    dx = 0.1;
    dy = dx;
    
    int = round(2*lx/dx)+1.0;
    l = round(lr/dh);
    k = 6;
    
    % The loading of the smoothened potentials for computing the
    % energy.    
    data = load(['./Smoothened_Potentials_New/3smoothened_logpotential_phi' num2str(phi(k)) '_stress' num2str(stress(pts)) '_check.txt']);
   
    % The initialization of x and y range of the potential and the potential itself.      
    xrange = [-lx1:dx:lx1];
    yrange = [-ly1:dy:ly1];
    V_xy = zeros(length(xrange),length(yrange));
    
    %%% While reading the potentials, one row is filled at a time since potential files named "...check.txt" is written first with 
        %%% fixed x and varying y.
    for bx = 1:1:length(xrange)
        V_xy(bx,[1:1:length(yrange)]') = data((bx-1.000)*length(xrange)+[1:1:length(yrange)]',3);
    end
    
    % Physical Dimensions of the hard core potential  
    Hard_Repulsion = 10^20;
    b = 3.0;
    sig = 6.0;
    
     % Initialization of the seed
    display('seed')
    seed = rnd_numbs(pts1);
    
    % The physical variables of the weights for the box fluctuations    
    fp =  0.00065;
    mu1 = 0.0;
    mu2 = 0.0;
    mu3 = 0.0;
    alp = 0.0;
    
    % Initial box parameters
    aspect = mu(pts1);
    alpha = 15.000;
    
    % Initial Box shape dimensions
    F_xx = -alpha;
    F_xy = alpha/aspect;
    F_yx = -alpha/aspect;
    F_yy = alpha;
    
    Area = abs((F_xx*F_yy)-(F_xy*F_yx));
    
    % Number of particles
    n = 3.0;
    N = NN(n);
    
    display(N)
    
    % The vertices of the initial box (x,y) coordinates
    x1 = 0;
    y1 = 0;
    x2 = 1*F_xx;
    y2 = 1*F_xy;
    x3 = 1*F_xx+1*F_yx;
    y3 = 1*F_xy+1*F_yy;
    x4 = 1*F_yx;
    y4 = 1*F_yy;
    
    xv = [x1 x2 x3 x4 x1]';
    yv = [y1 y2 y3 y4 y1]';
    
    % Initialization of the positions of particles inside the box     
    pos = zeros(N,2);
    pos_new2 = zeros(2,N);    

    % The positions of N particles is assigned here 
    pos_new2(1,:) = rand(1,N);
    pos_new2(2,:) = rand(1,N);
    
    % The affine transformation of positions of particles into a
    % square box
    h = [F_xx F_yx; F_xy F_yy]*pos_new2;
    
    pos(:,1) = h(1,:)';
    pos(:,2) = h(2,:)';
    
    % The number of metropolis moves at each bin of mu 
    MCS = 10^6;
   
    % Initialization of Required Variables
    Mac_Var = [];
    pos1 = [];
    Energy = 0;
    E = zeros(MCS,1);
    Acc = zeros(MCS,1);
    Acc_Global = zeros(MCS,1);
    dE_Global = zeros(MCS,1);
    rho0 = zeros(MCS,1);
    
    % Initialization of the correlation functions
    Cor1 = zeros(l,1);
    Cor2D = zeros(int,int);
    
    G1 = zeros(l,1);
    G12D = zeros(int,int);
    
    % The coordinates (qx,qy) of periodic vectors in all 9 directions
    qx = [(-F_xx+[-1 0 1]*F_yx) (0+[-1 0 1]*F_yx)  (F_xx+[-1 0 1]*F_yx)];
    qy = [(-F_xy+[-1 0 1]*F_yy) (0+[-1 0 1]*F_yy)  (F_xy+[-1 0 1]*F_yy)];
    
    % The loop which finds the minimum distance of all pairwise
    % bonds in one direction with the corresponding bonds in the periodic copies and
    % calculates energy based on the minimum distance
    for particle1 = 1:1:(N-1)
        
        delx = pos((particle1+1):1:end,1)-(pos(particle1,1)*ones(length(pos((particle1+1):1:end,1)),1));
        dely = pos((particle1+1):1:end,2)-(pos(particle1,2)*ones(length(pos((particle1+1):1:end,2)),1));
        
        delx = delx+qx;
        dely = dely+qy;
        
        Dist = sqrt(delx.^2 + dely.^2);
        
        % Gets the smallest distance
        delx = sum(delx.*logical(Dist == min(Dist,[],2)),2);
        dely = sum(dely.*logical(Dist == min(Dist,[],2)),2);
        Dist = sum(Dist.*logical(Dist == min(Dist,[],2)),2);
        
        % Gets the grid (bx,by) from the pairwise distances which is used to compute the energy            
        bx = round((lx1+delx)/dx)+1.00;
        a = find(logical(bx<=length(V_xy(:,1)) & bx>=1)==0);
        bx(a) = length(V_xy(:,1));
        
        by = round((ly1+dely)/dy)+1.00;
        a = find(logical(by<=length(V_xy(1,:)) & by>=1)==0);
        by(a) = length(V_xy(1,:));
        
        %%% While calculating the energy, one has to ensure, the matrix for potential is labelled by a single linear index by incrementing by one 
        %%% while going down a column.
        Energy = Energy+sum((V_xy(((by-1.000)*length(V_xy(:,1)))+bx)+(exp((0.02./Dist).^2)-1)).*(1./(1+exp(b*(Dist-sig)))));
        
        % Identifying the grid for the radial correlation function  
        br = round(Dist/dh) + 1;
        a = find(logical(br<=length(G1(:,1)) & br>=1)==0);
        br(a) = length(G1(:,1));
        
        % The radial correlation function is updated
        G1(br,1) = G1(br,1)+(2.000./(Dist));
        
        % Identifying the grid for the 2D correlation function  
        Dx = round((delx+lx)/(dx))+1.000;
        a = find(logical(Dx<=length(G12D(:,1)) & Dx>=1)==0);
        Dx(a) = length(G12D(:,1));
        
        Dy = round((dely+ly)/(dy))+1.000;
        a = find(logical(Dy<=length(G12D(:,1)) & Dy>=1)==0);
        Dy(a) = length(G12D(:,1));
        
        G12D((Dy-1.000)*int+Dx) = G12D((Dy-1.000)*int+Dx) + 1.000;
        
        Dx = round((-delx+lx)/(dx))+1.000;
        a = find(logical(Dx<=length(G12D(:,1)) & Dx>=1)==0);
        Dx(a) = length(G12D(:,1));
        
        Dy = round((-dely+ly)/(dy))+1.000;
        a = find(logical(Dy<=length(G12D(:,1)) & Dy>=1)==0);
        Dy(a) = length(G12D(:,1));
        
        G12D((Dy-1.000)*int+Dx) = G12D((Dy-1.000)*int+Dx) + 1.000;
        
    end
    
    % Double the energy for getting energy from all pairwise bonds in both directions
    Energy1 = (2*Energy);
   
    % Inputting the temperature for each parallel processing
    beta = b1(pts1);
    
    ang = 0.0;
    ang1 = 0.0;
    
    % The timer for the code starts from here.
    tic
    
    % Initialization of Required Variables
    ind = 0;
    step = 1;
    Accept = 0;
    Accept_Global = 0;
    i = 0;
    
    % The co-ordinates along x and y directions
    x = -lx:dx:lx;
    y = -ly:dy:ly;
    
    [X,Y] = meshgrid(x,y);

    % The initialization of variables that is needed to be stored.
    pos_new = zeros(N,2);
    dF_xy = 0;
    dF_yx = 0;
    dA = 0;
    posnew = [];
    posnew1 = [];   
    array1 = [];
    array2 = [];
    ix = 1; 
    iy = 1;

    % The loop begins which does the important sampling for N*(Monte Carlo
    % Steps)
    while(i<=N*MCS)
        
        i=i+1;
        
        if(mod(i,N)==0)
            
            stress(pts)
            N;
            fp;
            ind
            Energy1;
            F_xy
            -fp*N*Area;
            Energy1+fp*N*Area;
            Accept/i;
            Accept_Global/ind;
            
        end
        
        pos_new = pos;
        % The particle that is needed to be updated is chosen here
        particle = randi(N);
        
        pos_new(particle,1) = pos(particle,1);
        pos_new(particle,2) = pos(particle,2);
       
        % The positions of particles scaled back to the parallelogram with
        % given periodic vectors
        h = ([F_yy -F_yx; -F_xy F_xx]*[pos_new(particle,1) pos_new(particle,2)]')/Area;
          
        % We are producing a gaussian random variable for displacing the particles. 
        s1 = 1/sqrt(beta*N);
        r1 = s1*sqrt(-2*log(rand(1)));
        %r1 = s1*sqrt(rand(1));
        % The production of uniform angular distribution 
        theta = 2*pi*rand(1);
	
        array1 = [array1; r1 theta];

        % Gaussian random variables which are produced  for displacement
        new_pos_x = r1*cos(theta);
        new_pos_y = r1*sin(theta);
        
        % The displacement is added to existing positions
        h(1,1) = h(1,1)+new_pos_x;
        h(2,1) = h(2,1)+new_pos_y;
        
        % The particles are made to put into the square box by using a mod function.
        hx = mod(h(1,1),1);
        hy = mod(h(2,1),1);

        
        ix = ix*sign(hx);
        iy = iy*sign(hy);

	if(ix==-1)
		display('ix less than 1')
		hx
		hy
	end

	if(iy==-1)
		display('iy less than 1')
		hx
		hy
    end

        % The new positions of the particle are transformed to the
        % parallelogram box and relabelled
        h = [F_xx F_yx; F_xy F_yy]*[hx hy]';
       
        pos_new(particle,1) = h(1,1);
        pos_new(particle,2) = h(2,1);
        
        % The correlation function is defined here.
        G1_old = zeros(l,1);
        G12D_old = zeros(int,int);
        G1_new = zeros(l,1);
        G12D_new = zeros(int,int);
       
        % Gets the smallest distance for the particles in the new
        % configuration
        delx = pos_new(:,1)-(pos_new(particle,1)*ones(length(pos_new(:,1)),1));
        dely = pos_new(:,2)-(pos_new(particle,2)*ones(length(pos_new(:,2)),1));
        delx(particle) = [];
        dely(particle) = [];
        
        delx = delx+qx;
        dely = dely+qy;
        
        Dist = sqrt(delx.^2 + dely.^2);
        
        delx = sum(delx.*logical(Dist == min(Dist,[],2)),2);
        dely = sum(dely.*logical(Dist == min(Dist,[],2)),2);
        Dist = sum(Dist.*logical(Dist == min(Dist,[],2)),2);
        
        % Gets the grid (bx,by) from the pairwise distances which is used to compute the energy  
        bx = round((lx1+delx)/dx)+1.00;
        a = find(logical(bx<=length(V_xy(:,1)) & bx>=1)==0);
        bx(a) = length(V_xy(:,1));
        
        by = round((ly1+dely)/dy)+1.00;
        a = find(logical(by<=length(V_xy(1,:)) & by>=1)==0);
        by(a) = length(V_xy(1,:));
        
        % The energy of the new configuration.
        %%% While calculating the energy, one has to ensure, the matrix for potential is labelled by a single linear index by incrementing by one 
                %%% while going down a column.
        Energy_new = sum((V_xy(((by-1.000)*length(V_xy(:,1)))+bx)+(exp((0.02./Dist).^2)-1)).*(1./(1+exp(b*(Dist-sig)))));
        
        % Gets the grid for radial correlation function for new
        % configuration
        br = round(Dist/dh)+1.00;
        a = find(logical(br<=length(G1_new(:,1)) & br>=1)==0);
        br(a) = length(G1_new(:,1));
        
        % Computes the radial correlation function of new configuration
        G1_new(br,1) = G1_new(br,1)+(2.000./(Dist));
        
        % Computes the grid for 2D Correlation function
        Dx = round((delx+lx)/(dx))+1.000;
        a = find(logical(Dx<=length(G12D_new(:,1)) & Dx>=1)==0);
        Dx(a) = length(G12D_new(:,1));
        
        Dy = round((dely+ly)/(dy))+1.000;
        a = find(logical(Dy<=length(G12D_new(1,:)) & Dy>=1)==0);
        Dy(a) = length(G12D_new(:,1));
        
        % Computes the 2D correlation function of new configuration
        G12D_new((Dy-1.000)*int+Dx) = G12D_new((Dy-1.000)*int+Dx) + 1.000;
        
        
        Dx = round((-delx+lx)/(dx))+1.000;
        a = find(logical(Dx<=length(G12D_new(:,1)) & Dx>=1)==0);
        Dx(a) = length(G12D_new(:,1));
        
        Dy = round((-dely+ly)/(dy))+1.000;
        a = find(logical(Dy<=length(G12D_new(1,:)) & Dy>=1)==0);
        Dy(a) = length(G12D_new(:,1));
        
        G12D_new(((Dy-1.000)*int)+Dx) = G12D_new(((Dy-1.000)*int)+Dx) + 1.000;
        
         % Gets the smallest distance for the particles in the old
        % configuration
        delx = pos(:,1)-(pos(particle,1)*ones(length(pos(:,1)),1));
        dely = pos(:,2)-(pos(particle,2)*ones(length(pos(:,2)),1));
        delx(particle) = [];
        dely(particle) = [];
        
        delx = delx+qx;
        dely = dely+qy;
        
        Dist = sqrt(delx.^2 + dely.^2);
        
        delx = sum(delx.*logical(Dist == min(Dist,[],2)),2);
        dely = sum(dely.*logical(Dist == min(Dist,[],2)),2);
        Dist = sum(Dist.*logical(Dist == min(Dist,[],2)),2);
        
         % Gets the grid (bx,by) from the pairwise distances which is used to compute the energy  
        bx = round((lx1+delx)/dx)+1.000;
        a = find(logical(bx<=length(V_xy(:,1)) & bx>=1)==0);
        bx(a) = length(V_xy(:,1));
        
        by = round((ly1+dely)/dy)+1.000;
        a = find(logical(by<=length(V_xy(1,:)) & by>=1)==0);
        by(a) = length(V_xy(1,:));
        
         % The energy of the old configuration.
        %%% While calculating the energy, one has to ensure, the matrix for potential is labelled by a single linear index by incrementing by one 
                %%% while going down a column.
        Energy_old = sum((V_xy(((by-1.000)*length(V_xy(:,1)))+bx)+(exp((0.02./Dist).^2)-1)).*(1./(1+exp(b*(Dist-sig)))));
        
        % Identifying the grid for the radial correlation function  
         
        br = round(Dist/dh)+1.00;
        a = find(logical(br<=length(G1_old(:,1)) & br>=1)==0);
        br(a) = length(G1_old(:,1));
        
        % Computes the radial correlation function of old configuration
        G1_old(br,1) = G1_old(br,1)+(2.0000./(Dist));
        
        Dx = round((delx+lx)/(dx))+1.000;
        a = find(logical(Dx<=length(G12D_old(:,1)) & Dx>=1)==0);
        Dx(a) = length(G12D_old(:,1));
        
        Dy = round((dely+ly)/(dy))+1.000;
        a = find(logical(Dy<=length(G12D_old(1,:)) & Dy>=1)==0);
        Dy(a) = length(G12D_old(:,1));
        
         % Computes the 2D correlation function of old configuration
        G12D_old((Dy-1.000)*int+Dx) = G12D_old((Dy-1.000)*int+Dx) + 1.000;
        
        Dx = round((-delx+lx)/(dx))+1.000;
        a = find(logical(Dx<=length(G12D_old(:,1)) & Dx>=1)==0);
        Dx(a) = length(G12D_old(:,1));
        
        Dy = round((-dely+ly)/(dy))+1.000;
        a = find(logical(Dy<=length(G12D_old(1,:)) & Dy>=1)==0);
        Dy(a) = length(G12D_old(:,1));
        
        G12D_old(((Dy-1.000)*int)+Dx) = G12D_old(((Dy-1.000)*int)+Dx) + 1.000;
        
        % Modifies the total energy of the system.
        Energy3 = Energy1-(2*Energy_old)+(2*Energy_new);
        
        % This loop is executed based on boltzmann weight using difference in
        % energy.
        if(rand(1) < exp(-beta*(Energy3-Energy1)))
            % The positions, energy, radial, 2D Correlation Function, Accept Ratio is modified here 
            pos = pos_new;
            Energy1 = Energy3;
            Accept = Accept+1;
            G1 = G1-G1_old+G1_new;
            G12D = G12D-G12D_old+G12D_new;
        end
        
        % The quantities are evaluated after every N steps which is one
        % Monte Carlo step.
        if(mod(i,N)==0)
            
            ind = ind+1;
            
            % The Energy, 2D Correlation Function, Radial Correlation Fn is
            % updated here.
            E(ind,1) = Energy1/N;
            % conf(ind,1) = {ind};
            % conf(ind,2) = {[xv yv]};
            % conf(ind,3) = {F_xy};
            % conf(ind,4) = {F_yx};
            % conf(ind,5) = {pos};
            Cor2D = Cor2D+(G12D./(N*dx*dy*(N/Area)));
            Cor1 = Cor1+(G1./(N*2*pi*dh*(N/Area)));
            
     
            
            % The macro variables are stored in this array.
           Mac_Var = [Mac_Var; (Energy1/N) F_xx F_xy abs(F_xx/F_xy)];
            
           % The quantities are evaluated after every 100
            % Monte Carlo step.       
           if(mod(ind,100)==0)
                 pos1 = [pos1; {pos}];
           end

            % The random number is decided based on which the periodic
            % vector is chosen that needs to be modified.
            rndm = randi(2);
            % The old area is stored in this variable.
            A_old = Area ;
            
            if(rndm==1)
                
                done = 1;
                
                % F_xy is updated here and the particles are transfered
                % here into the new box.
                while(done<=1)
                      
                     % The increment size for periodic vector is chosen here.
                     dA = 2*(rand()-0.5)*(dstep(pts1));
                     dF_xy = sqrt((A_old+dA)+(F_xx)^2)-sqrt(A_old+(F_xx)^2) ;	
                   
                    % An if condition to put a threshold on how small the
                    % box can get
                    if(abs(F_xy+dF_xy) >= abs(F_xx) && abs(F_yx-dF_xy) >= abs(F_yy))
                        
                        % The positions of particles after updation are scaled into a square box
                        % of unit area
                        posnew1 = ([F_yy -F_yx; -F_xy F_xx]*[pos(:,1) pos(:,2)]')/Area;
                        %A = [1.0 0.0 ; -(F_yy*dF_xy/((F_yx*F_xy)-(F_xx*F_yy))) 1.0+(F_yx*dF_xy/((F_yx*F_xy)-(F_xx*F_yy)))];
                        
                        % The periodic vectors are updated here
                        F_xy = F_xy+dF_xy;
                        
                        %posnew1 = A*pos';
                        %posnew1 = posnew1';
                        
                        dF_yx = -dF_xy;
                        
                        %B = [1.0+(F_xy*dF_yx/((F_xy*F_yx)-(F_yy*F_xx))) -(F_xx*dF_yx/((F_xy*F_yx)-(F_yy*F_xx))); 0.0 1.0];
                        
                        %posnew = B*posnew1';
                        %posnew = posnew';
                        
                        F_yx = F_yx+dF_yx;
                        
                        % The positions of particles scaled earlier in this loop 
                        % to square box is transformed back to
                        % parallelogram with updated periodic vectors.
                        posnew = [F_xx F_yx; F_xy F_yy]*posnew1;
                        posnew = posnew';
                        
                        done = done+1;
                        
                    end
                end
                else
                
                done = 1;
                
                % F_yx is updated here and the particles are transfered
                % here into the new box.
                while(done<=1)
                    
                    % The positions of particles after updation are scaled into a square box
                        % of unit area
                     dA = 2*(rand()-0.5)*(dstep(pts1));
                     dF_yx = -(sqrt((A_old+dA)+(F_xx)^2)-sqrt(A_old+(F_xx)^2)) ;	
                    
                      % An if condition to put a threshold on how small the
                    % box can get
                    if(abs(F_yx+dF_yx) >= abs(F_yy) && abs(F_xy-dF_yx) >= abs(F_xx))
                        
                        % The positions of particles after updation are scaled into a square box
                        % of unit area
                        posnew1 = ([F_yy -F_yx; -F_xy F_xx]*[pos(:,1) pos(:,2)]')/Area;
                        %B = [1.0+(F_xy*dF_yx/((F_xy*F_yx)-(F_yy*F_xx))) -(F_xx*dF_yx/((F_xy*F_yx)-(F_yy*F_xx))); 0.0 1.0];
                        
                        %posnew1 = B*pos';
                        %posnew1 = posnew1';
                        
                       % The periodic vectors are updated here
                        F_yx = F_yx+dF_yx;
                        
                        dF_xy = -dF_yx;
                        %A = [1.0 0.0 ; -(F_yy*dF_xy/((F_yx*F_xy)-(F_xx*F_yy))) 1.0+(F_yx*dF_xy/((F_yx*F_xy)-(F_xx*F_yy)))];
                        
                        %posnew = A*posnew1';
                        %posnew = posnew';
                        
                        F_xy = F_xy+dF_xy;
                        
                        % The positions of particles scaled earlier in this loop 
                        % to square box is transformed back to
                        % parallelogram with updated periodic vectors.
                        posnew = [F_xx F_yx; F_xy F_yy]*posnew1;
                        posnew = posnew';
                        
                        done = done+1;
                    end
                end
            end
            
            % The array stores the incremented area (dA) after each global
            % move.
            array2 = [array2; dA];

            % The coordinates (qx,qy) of periodic vectors in all 9 directions 
            qx = [(-F_xx+[-1 0 1]*F_yx) (0+[-1 0 1]*F_yx)  (F_xx+[-1 0 1]*F_yx)];
            qy = [(-F_xy+[-1 0 1]*F_yy) (0+[-1 0 1]*F_yy)  (F_xy+[-1 0 1]*F_yy)];
            
           % Initialization of energy
            Energy = 0;
            
            % Initialization of correlation functions that is updated
            % after one global moves.
            G1_Global = zeros(l,1);
            G12D_Global = zeros(int,int);
            
             % The loop which finds the minimum distance of all pairwise
            % bonds in one direction with the corresponding bonds in the periodic copies and
            % calculates energy based on the minimum distance
            for particle1=1:1:(N-1)
                
                delx = posnew((particle1+1):1:end,1)-(posnew(particle1,1)*ones(length(posnew((particle1+1):1:end,1)),1));
                dely = posnew((particle1+1):1:end,2)-(posnew(particle1,2)*ones(length(posnew((particle1+1):1:end,2)),1));
                
                delx = delx+qx;
                dely = dely+qy;
                
                Dist = sqrt(delx.^2 + dely.^2);
                
                % Gets the smallest distance
                delx = sum(delx.*logical(Dist == min(Dist,[],2)),2);
                dely = sum(dely.*logical(Dist == min(Dist,[],2)),2);
                Dist = sum(Dist.*logical(Dist == min(Dist,[],2)),2);
                
                % Gets the grid (bx,by) from the pairwise distances which is used to compute the energy  
                bx = round((lx1+delx)/dx)+1.00;
                a = find(logical(bx<=length(V_xy(:,1)) & bx>=1)==0);
                bx(a) = length(V_xy(:,1));
                
                by = round((ly1+dely)/dy)+1.00;
                a = find(logical(by<=length(V_xy(1,:)) & by>=1)==0);
                by(a) = length(V_xy(1,:));
                
                %%% While calculating the energy, one has to ensure, the matrix for potential is labelled by a single linear index by incrementing by one 
                %%% while going down a column.

                Energy = Energy + sum((V_xy(((by-1.000)*length(V_xy(:,1)))+bx)+(exp((0.02./Dist).^2)-1)).*(1./(1+exp(b*(Dist-sig)))));
                
                % Gets the grid for computing radial pair correlation
                % functions
                br = round(Dist/dh)+1.0;
                a = find(logical(br<=length(G1_Global(:,1)) & br>=1)==0);
                br(a) = length(G1_Global(:,1));
                
                % The radial correlation functions gets updated after one
                % global move
                G1_Global(br,1) = G1_Global(br,1)+(2.000./(Dist));
                
                % Computes the grid for 2D Correlation function
                Dx = round((delx+lx)/(dx))+1.000;
                a = find(logical(Dx<=length(G12D_Global(1,:)) & Dx>=1)==0);
                Dx(a) = length(G12D_Global(1,:));
                
                Dy = round((dely+ly)/(dy))+1.000;
                a = find(logical(Dy<=length(G12D_Global(1,:)) & Dy>=1)==0);
                Dy(a) = length(G12D_Global(1,:));
                
                % The 2D correlation functions gets updated after one
                % global move
                G12D_Global((Dy-1.000)*int+Dx) = G12D_Global((Dy-1.000)*int+Dx) + 1.000;
                
                Dx = round((-delx+lx)/(dx))+1.000;
                a = find(logical(Dx<=length(G12D_old(:,1)) & Dx>=1)==0);
                Dx(a) = length(G12D_old(:,1));
                
                Dy = round((-dely+ly)/(dy))+1.000;
                a = find(logical(Dy<=length(G12D_old(1,:)) & Dy>=1)==0);
                Dy(a) = length(G12D_old(1,:));
                
                G12D_Global((Dy-1.000)*int+Dx) = G12D_Global((Dy-1.000)*int+Dx) + 1.000;
            end
            
            % The new area is stored here
            A_new = abs((F_xx*F_yy)-(F_yx*F_xy));
            % dE_Global(ind,2) = (A_new-A_old)/A_new;
            % The energy for considering bonds from both direction is therefore doubled here 
            Energy3 = (2*Energy);
            
           % display('Prob')
            c = rand(1);
            %display('Weight')
            %exp(-b1*(Energy3-Energy1)-(b1*fp*N*(A_new-A_old))+(b1*N*log(A_new/A_old)));
            
            %dE_Global(ind,1) = exp(-b1*(Energy3-Energy1));
            %dE_Global(ind,3) = Energy3-Energy1;
            
            % This loop is executed based on boltzmann weight using difference in
        % energy, area, and log of area.
            if(c < exp(-(beta*(Energy3-Energy1))-(beta*fp*N*(A_new-A_old))+(N*(log(A_new)-log(A_old)))))
                % The positions, energy, radial, 2D Correlation Function, Energy, Accept Ratio is modified here 
                pos = posnew;
                Energy1 = Energy3;
                Accept = Accept+1;
                Accept_Global = Accept_Global+1;
                
                G1 = G1_Global;
                G12D = G12D_Global;
                
                x1 = 0.0;
                y1 = 0.0;
                x2 = 1.0*F_xx;
                y2 = 1.0*F_xy;
                x3 = 1.0*F_xx+1.0*F_yx;
                y3 = 1.0*F_xy+1.0*F_yy;
                x4 = 1.0*F_yx;
                y4 = 1.0*F_yy;
                
                xv = [x1 x2 x3 x4 x1]';
                yv = [y1 y2 y3 y4 y1]';
                
                Area = A_new;
            else
                % We go back to the older configurations
                F_xy = F_xy-dF_xy;
                F_yx = F_yx-dF_yx;
                Area = A_old;
                qx = [(-F_xx+[-1 0 1]*F_yx) (0+[-1 0 1]*F_yx) (F_xx+[-1 0 1]*F_yx)];
                qy = [(-F_xy+[-1 0 1]*F_yy) (0+[-1 0 1]*F_yy) (F_xy+[-1 0 1]*F_yy)];
            end
            
            % The accept ratio for local anf global moves is scaled
            Acc(ind,1) = Accept/(ind*(N+1));
            %rho0(ind,1) = N/abs((F_xx*F_yy)-(F_xy*F_yx)) ;
            Acc_Global(ind,1) = Accept_Global/ind;
            
            % After pre-decided number of Monte Carlo steps relevant quanties are saved into file
            if(ind == set(step))
                Var_E = zeros(ind,1);
                E_avg=0;
                E2_avg =0;
                
                for k1=1:1:ind
                    E_avg = E_avg+(E(k1,1));
                    E2_avg = E2_avg+(E(k1,1)*E(k1,1)/(1*1));
                    Var_E(k1,1) = sqrt(((E2_avg/k1)-((E_avg/k1)^2)));
                end
                
                %  file = ['Conf_N64_NpT_ensemble_fix_Gxy_' num2str(s)];
                %         parsave(file,E(1:ind,:),Acc(1:ind,:),Var_E(1:ind,:)...
                % ,conf(1:ind,:),rho0(1:ind,:),Acc_Global(1:ind,:),dE_Global(1:ind,:));
                
                
                Ratio = [Var_E(1:ind,:) Acc(1:ind,:) Acc_Global(1:ind,:)];
                
                % Creating a file for storing required values 
                file = ['./NewResults/Cor_ICV_Gaussian1_sig6_sample' num2str(pts1)  '_beta' num2str(b1(pts1)) '_smallstep' num2str(dstep(pts1)) '_mu0.33_phi' num2str(phi(k)) '_stress' num2str(stress(pts)) '.mat'];
                parsave1(file,Cor1./ind,Cor2D./ind, pos1, Mac_Var, Ratio, array1, array2);
                
                step = step+1;
            end
        end
    end
    % The timer is switched off.
    toc
end


% Saving the table of relevant variables in the created file. 
function parsave1(file,C1,C2D,pos1,Mac_Var,Ratio,array1,array2)
save(file,'C1','C2D','pos1','Mac_Var','Ratio','array1','array2');
end

