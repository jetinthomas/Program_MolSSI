% The initialization and number of cores needed in the execution of the code 
parpool(6);
% The number of particles that is in the simulation box 
NN = [100 128 256 512 768 1024 1280 1408 1536 3000];
% The strength of the hard-core potential multiplied to the smoothened
% potential to remove the intense clustering at very small scales.
vo = [1 2 3 4 5 10 20 50 100 200 500 1000 2000 5000 10000];
% for loop for going through different cores 
parfor n = 1:1:6
    
    % The different stresses in the phase diagram for DST suspensions 
    stress = [0.1 0.5 1.0 2.0 5.0 10.0 20.0 100.0];
    
    % The packing fractions in the phase diagram for DST suspensions
    phi = [0.76 0.77 0.78 0.785 0.79 0.8 0.82];
    
    display('N')
    display(n);

    % Initializing the strength of the Hard Core potential
    display('ep')
    ep = vo(15);
    display(ep);
    
    % Use of different random numbers for initializing the seed for different parallel cores. 
    rng('shuffle','twister');
    rnd_numbs = (10^9).*rand(8,1)
    
   % for loop for going through different cores 
    for pts = 8:1:8
        
       % set = [10 100 1000 10000 100000 200000 500000 800000 1000000 2000000 4000000 6000000 8000000 10000000 20000000 50000000 80000000 100000000];
        
       % Useful dimensions required for the program
        lr = 50.0;
        lx = 50.0;
        ly = 50.0;
        lx1 = 10.0;
        ly1 = 10.0;
        
        dhx = 0.1;
        dhy = 0.1;
        
        dh = 0.2;
        
        dx = 0.1;
        dy = dx;
        
        int = round(2*lx/dx);
        l = round(lr/dh);
        
        k=6;
        
        % The loading of the smoothened potentials for computing the
        % energy.
        data = load(['./Smoothened_Potentials_New/3smoothened_logpotential_phi' num2str(phi(k)) '_stress' num2str(stress(pts)) '_check.txt']);
       
        % The initialization of x and y range of the potential and the potential itself. 
        xrange = -lx1:dx:lx1;
        yrange = -ly1:dy:ly1;
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
        seed = rnd_numbs(n);

        % The physical variables of the weights for the box fluctuations
        fp =  0.0;
        mu1 = 0.0;
        mu2 = 0.0;
        mu3 = 0.0;
        alp = 1.0;
        
        % Initial box parameters
        aspect = 1.0000/3.0000;
        alpha = 15.000;
        
        % Box Area (1/rho) Increments 
        dr = 0.000005;
        % Starting Box Area (1/rho)
        rho = 0.000005;
        
        % Initial Box shape dimensions
        F_xx = -alpha;
        F_xy = sqrt((1/rho)+225);
        F_yx = -F_xy;
        F_yy = alpha;
        
        % Number of particles
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
        pos_new = zeros(2,N);
        
        % The number of metropolis moves at each bin of mu 
        MCS = 400;
       
        % Initialization of the table which stores all the required values. 
        E = zeros(500,8);
        E(:,6:7) = 100000*ones(500,2);
        
        Energy1 = 0;
        
        % The timer for the code starts from here.
        tic
        
        % The initialization of relevant flags
        ind = 1;
        ind1 = 0;
        ind2 = 0;
        step = 1;
        i = 1;
        
        aa = 1.0;
        
        % The loop that goes from the smallest rho to largest rho and
        % backwards several times
        while(ind2 <= MCS*820)
            
            % The positions of N particles is assigned here in a square
            % box.
            pos_new(1,:) = rand(1,N);
            pos_new(2,:) = rand(1,N);
            
            % The affine transformation of positions of particles into a
            % parallelogram whose dimensions are given by periodic vectors 
            h = [F_xx F_yx; F_xy F_yy]*pos_new;
            
            pos(:,1) = h(1,:)';
            pos(:,2) = h(2,:)';
            
            % Initialization of energy
            Energy = 0;
            
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
                bx = round((lx1+delx)/dx)+1;
                a = find(logical(bx<=length(V_xy(:,1)) & bx>=1)==0);
                bx(a) = length(V_xy(:,1));
                
                by = round((ly1+dely)/dy)+1;
                a = find(logical(by<=length(V_xy(1,:)) & by>=1)==0);
                by(a) = length(V_xy(:,1));
                
                %%% While calculating the energy, one has to ensure, the matrix for potential is labelled by a single linear index by incrementing by one 
                %%% while going down a column.
 
                Energy = Energy + sum((ep*(V_xy(((by-1.000)*length(V_xy(:,1)))+bx))+(exp((0.02./Dist).^2)-1.0)).*(1./(1+exp(b*(Dist-sig)))));
                
            end
            
            % Double the energy for getting energy from all pairwise bonds in both directions
            Energy1 = (2*Energy);
            
            % Area of the force tile box given the periodic vectors
            Area = abs((F_xx*F_yy)-(F_xy*F_yx));
            
            % Put an threshold on energy
            if(Energy1 < 100000)
                
                % Incrementing the values at the corresponding rho bin
                E(round(rho/dr),1) = E(round(rho/dr),1)+ 1.0;
                E(round(rho/dr),2) = E(round(rho/dr),2)+ (N*log(Area));
                E(round(rho/dr),3) = E(round(rho/dr),3)+ exp(-Energy1);
                E(round(rho/dr),4) = E(round(rho/dr),4)+ Energy1;
                E(round(rho/dr),5) = rho;
                
                if(Energy1 < E(round(rho/dr),6))
                    E(round(rho/dr),6) = Energy1;
                    E(round(rho/dr),7) = Energy1;
                end
                
                ind1 = ind1+1.0;
                
            else
                
                E(round(rho/dr),8) =  E(round(rho/dr),8)+1.0;
                
            end
            
            % Printing certain values after each monte carlo steps
            if(ind1 == MCS)
                
                display('After MCS');
                phi(k)
                stress(pts)
                N
                i
                ind1
                ind2
                % F_xy
            	ep
                sprintf('%1.10f',rho)
           
            end
            
            % The Area update after each monte carlo steps
            if(ind1 == MCS)
                
                done = 1;
                
                while (done <=1)
                    
                    % The minimum rho below which you cant go.
                    if(rho <= 0.0000055)
                        aa = 1;
                    end
                    
                    % The maximum rho above which you cant go
                    if(rho >= 0.002000)
                        aa = -1;
                    end
                    
                    % The increment in rho and periodic vector
                    drho = aa*0.000005;
                    dF_xy = sqrt((1/(rho+drho))+225)-sqrt((1/rho)+225);
                    
                    % An if condition to put a threshold on how small the
                    % box can get
                    if(abs(F_xy+dF_xy) >= abs(F_xx) && abs(F_yx-dF_xy) >= abs(F_yy))
                        % The box dimensions are incremented here
                        F_xy = F_xy+dF_xy;
                        dF_yx = -dF_xy;
                        F_yx = F_yx+dF_yx;
                        rho = rho+drho;
                        done = done+1;
                    end
                end
                
                % A Monte Carlo number of steps is added to global variable
                % 
                ind2 = ind2+ind1;
                ind1 = 0;
                
                EE = [E(:,6) E(:,8) E(:,1)];
                
                % Creating a file for storing required values 
                file = ['./pf80_s100_Analysis/E_min_' num2str(N) '_smallI_randbox_vo' num2str(ep) '_phi' num2str(phi(k)) '_stress' num2str(stress(pts)) '.mat'];
                parsave(file, seed, EE);
                
            end
            i=i+1;
        end
        
        % The timer is switched off.
        toc
        
    end
end

% Saving the table of relevant variables along with the seed. 
 function parsave(file,seed,EE)
 save(file,'seed','EE');
 end
