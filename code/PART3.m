  % PART 3: PUTTING IT ALL TOGETHER
    clearvars
    particles = 10;
    
    global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons

    %% SETTING UP THE CANVAS ----------------------------------------------
    
    % The nominal size of the region is 200nm × 100nm
    region_x = 200e-9;
    region_y = 100e-9;
    % Pick time step size which should be smaller than 1/100 the region
    step_size = 1e-9;
    timestep = 1000;
    % Electron Concentration is 10e15 cm^-2
    elecC = 10e15;
    % Assume T = 300K
    T = 300;
    % Voltage applied across x dimension of semiconductor
    V = 0.1;
    % Calculate thermal velocity
    v_th = sqrt(2*C.kb*T/C.m_n);
    % Change in velocity for each timestep
    v_change = step_size/v_th;
    % Mean time between collisions is 0.2e-12
    MT_C = 0.2e-12;
    % Conductivity values
    sigma_o = 1;        % outside boxes 
    sigma_in = 1e-2;    % inside boxes
    
    
    %% RANDOM PARTICLE INITIALIZATION -------------------------------------
    
    % Creates 2 rows of random numbers for X coordinates of each particle.
    % The two rows are used for position and angle.
    % Y coordinates only requires position.
    X = rand(2,particles);
    Y = rand(1,particles);
    % Set up the X and Y position coordinates for each particle
    % Multiply by region restrictions to ensure it's within the region
    X_position(1,:) = X(1,:)*region_x;
    Y_position(1,:) = Y(1,:)*region_y;
    
    % THESE PARTICLES CANNOT BE IN THE ADDED RECTANGLES
    % Need to check their current starting location, and if it's inside
    % the rectangles then we need a new random position.
    
    % Check if particle is inside X boundary
    check_X_left = X_position > 0.8e-7;
    check_X_right = X_position < 1.2e-7;
    check_X = check_X_left & check_X_right;
    % Check which box it's in by checking Y coordinates
    check_top = Y_position > 0.6e-7;
    check_bottom = Y_position < 0.4e-7;
    % Find out the box
    box_top = check_top & check_X;
    box_bottom = check_bottom & check_X;
    IN_A_BOX = box_top | box_bottom;
    
    % Now we start randomizing the positions of particles inside box
    while(sum(IN_A_BOX) > 0)
        
        temp_x = rand(1,sum(IN_A_BOX));
        temp_y = rand(1,sum(IN_A_BOX));
        % Randomize position
        X_position(IN_A_BOX) = temp_x*region_x;
        Y_position(IN_A_BOX) = temp_y*region_y;
        %X_position(IN_A_BOX) = rand*region_x;
        %Y_position(IN_A_BOX) = rand*region_y;
        
        % Check again
        check_X_left = X_position > 0.8e-7;
        check_X_right = X_position < 1.2e-7;
        check_X = check_X_left & check_X_right;
        check_top = Y_position > 0.6e-7;
        check_bottom = Y_position < 0.4e-7;
        box_top = check_top & check_X;
        box_bottom = check_bottom & check_X;
        IN_A_BOX = box_top | box_bottom;
    end
    % Once we leave while loop, we have no particles in the boxes
    
    
    % Normal Distribution of Velocity
    v_th_N = v_th/sqrt(2);
    X_velocity = v_th_N*randn(1,particles)*v_change;
    Y_velocity = v_th_N*randn(1,particles)*v_change;
    
    % Set up scattering percent
    PSCAT = 1 - exp(-v_change/MT_C);
    mfp_vec = zeros(1,particles);
    
    % Set up the electric field
    E_FIELD = V/region_x;
    E_FORCE = E_FIELD*C.q_0;
    acceleration = E_FORCE/C.m_n;
    vel_fromA = acceleration*(v_change^2);
    
    CURRENT_STORAGE = zeros(1,timestep);
    
    % Set up box dimensions
    box_X = [0.8e-7 1.2e-7];
    box_Y_bottom = [0 0.4e-7];
    box_Y_top = [0.6e-7 region_y];
    
    % Map the electric field in 100 different bins
    bins = 100;
    x_efield = region_x/bins;
    y_efield = region_y/bins;
    % Make a meshgrid for the 100 bins from 0 to end of x and y regions
    [dimensions_X, dimensions_Y] = meshgrid(0:x_efield:region_x,0:y_efield:region_y);
    
    % TIME TO COMBINE THIS WITH THE FINITE DIFFERENCE STUFF!
    cond = zeros(bins+1); % we want to have all 100 bins that we just set up, and have a conductivity for all of them
    
    % Figure out what parts of the dimensions_X and dimensions_Y are inside
    % the boxes
    
    % Find all x dimensions inside the boxes
    box_X_logic = dimensions_X >= box_X(1) & dimensions_X <= box_X(2);
    % Find all y dimensions inside the boxes
    box_Y_logic = dimensions_Y >= box_Y_top(1) | dimensions_Y <= box_Y_bottom(2);
    
    box_logic_final = box_X_logic & box_Y_logic;
    
    % set up the conductivity using the logic assignments of the boxes
    % Everything inside box:
    cond(box_logic_final) = sigma_in;
    % Everything outside box:
    cond(~box_logic_final) = sigma_o;
    
    % G MATRIX SETUP TIME
    FIELD_VOLTAGE = 0.8; % whatever is going on is visible at 0.8 V so let's try that
    
    G = sparse(bins+1);
    V = zeros(bins+1,1);
    map = @(i,j) j + (i - 1)*bins;
    
    % X direction
    for i=1:bins
        % Y direction
        for j=1:bins
            n = map(i,j);
            nxm = map(i-1,j);
            nxp = map(i+1,j);
            nym = map(i,j-1);
            nyp = map(i,j+1);

            % nx = 0 (in this case, V = Vo)
            % we also don't have boxes on the left and right edges so don't
            % need to check anything :)
            if i == 1
                G(n,:) = 0;
                G(n,n) = cond(i,j);
                V(n) = FIELD_VOLTAGE;
            % nx = L (in this case, V = 0)
            elseif i == bins
                G(n,:) = 0;
                G(n,n) = cond(i,j);
                V(n) = 0;
            % Y edges of canvas
            elseif j == 1
                G(n,:) = 0;
                G(n,nxm) = (cond(i-1,j) + cond(i,j))/2;
                G(n,nxp) = (cond(i+1,j) + cond(i,j))/2;
                G(n,nyp) = (cond(i,j+1) + cond(i,j))/2;
                G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nyp));
            elseif j == bins
                G(n,:) = 0;
                G(n,nxm) = (cond(i-1,j) + cond(i,j))/2;
                G(n,nxp) = (cond(i+1,j) + cond(i,j))/2;
                G(n,nym) = (cond(i,j-1) + cond(i,j))/2;
                G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nym));

            else          
                G(n,:) = 0;
                G(n,nxm) = (cond(i-1,j) + cond(i,j))/2;
                G(n,nxp) = (cond(i+1,j) + cond(i,j))/2;
                G(n,nyp) = (cond(i,j+1) + cond(i,j))/2;
                G(n,nym) = (cond(i,j-1) + cond(i,j))/2;
                G(n,n) = -(G(n,nxm) + G(n,nxp) + G(n,nyp) + G(n,nym));  
            end

        end
    end

% Solve for F using G\V = F
F = G\V;
    
% Set up a surf plot
surfs_up = zeros(bins);
for i = 1:bins
    for j = 1:bins
        n = map(i,j);
        surfs_up(i,j) = F(n);
    end
end

[EX, EY] = gradient(surfs_up*10e6);

% Set up accelerations for use in plotting
X_FORCE = EX*C.q_0;
Y_FORCE = EY*C.q_0;
acceleration_X = X_FORCE/C.m_n;
acceleration_Y = Y_FORCE/C.m_n;
vel_X_fromA = acceleration_X*(v_change^2);
vel_Y_fromA = acceleration_Y*(v_change^2);

%% Begin plotting the boxes! ---------------------------------------------
    
    % Go through loop for specified number of time steps
    for p = 1:timestep
        
        % ACCELERATE THE PARTICLES WITH OUR E FIELD VELOCITIES!
        for L = 1:bins
            for W = 1:bins
                % Figure out if we're accelerating or not
                logic_accelX = X_position < L*x_efield & X_position > (L-1)*x_efield;
                logic_accelY = Y_position < W*y_efield & Y_position > (W-1)*y_efield;
                
                X_velocity(logic_accelX) = X_velocity(logic_accelX) + vel_X_fromA(L,W);
                Y_velocity(logic_accelY) = Y_velocity(logic_accelY) + vel_Y_fromA(L,W);
            end
        end
        
        % Scattering electrons
        scatter = rand(1,particles);
        DID_I_SCATTER = scatter < PSCAT;
        
        % If PSCAT is greater than scatter, the particle scatters.
        % In this case, all the 1's in logical array DID_I_SCATTER
        % will be scattered and need to have recalculated stuff.
        % Get new random velocity values for scattering purposes:
        X_velocity_rand = v_th_N*randn(1,particles);
        Y_velocity_rand = v_th_N*randn(1,particles);
       
        X_velocity(DID_I_SCATTER) = X_velocity_rand(DID_I_SCATTER)*v_change;
        Y_velocity(DID_I_SCATTER) = Y_velocity_rand(DID_I_SCATTER)*v_change;

        % Keep track of how many particles scattered
        % to figure out mean free path
        mfp_vec(~DID_I_SCATTER) = mfp_vec(~DID_I_SCATTER)+step_size;
        % Anything that scattered will be set to 0
        mfp_vec(DID_I_SCATTER) = 0;
        
        
        % BOX REFLECTIONS -------------------------------------------------
        
       % BOTTOM BOX TEST
       % Check if particle is inside X boundary
        check_X_left = (X_position + X_velocity) > 0.8e-7;
        check_X_right = (X_position + X_velocity) < 1.2e-7;
        check_X = check_X_left & check_X_right;
        % Check Y boundary
        check_bottom = (Y_position + Y_velocity) < 0.4e-7;
        % Check bottom box
        box_bottom = check_bottom & check_X;
        X_velocity(box_bottom) = -1*X_velocity(box_bottom);
        % Check if we need to multiply Y component by -1 as well
        check_X_left = (X_position + X_velocity) > 0.8e-7 + step_size;
        check_X_right = (X_position + X_velocity) < 1.2e-7 - step_size;
        check_X = check_X_left & check_X_right;
        check_Y_bottombox_top = Y_position < 0.4e-7 - step_size;
        Y_BOX_BOTTOM = check_X & check_Y_bottombox_top;
        Y_velocity(Y_BOX_BOTTOM) = -1*Y_velocity(Y_BOX_BOTTOM);
        
        % TOP BOX TEST
        % Check if particle is inside X boundary
        check_X_left = (X_position + X_velocity) > 0.8e-7;
        check_X_right = (X_position + X_velocity) < 1.2e-7;
        check_X = check_X_left & check_X_right;
        % Check Y boundary
        check_top = (Y_position + Y_velocity) > 0.6e-7;
        % Check bottom box
        box_top = check_top & check_X;
        X_velocity(box_top) = -1*X_velocity(box_top);
        % Check if we need to multiply Y component by -1 as well
        check_X_left = (X_position + X_velocity) > 0.8e-7 + step_size;
        check_X_right = (X_position + X_velocity) < 1.2e-7 - step_size;
        check_X = check_X_left & check_X_right;
        check_Y_topbox_bottom = Y_position > 0.6e-7 + step_size;
        Y_BOX_TOP = check_X & check_Y_topbox_bottom;
        Y_velocity(Y_BOX_TOP) = -1*Y_velocity(Y_BOX_TOP);
        
        % Begin moving the electrons by using initial position
        % and adding the velocity to it
        % Use logical array to make stuff easier...
        
        % Keep all the particles within the canvas.
        % If their position + velocity goes over the right X boundary,
        % We will subtract region_x to make it jump to the opposite edge.
        X_right = X_position + X_velocity > region_x;
        X_position(X_right) = X_position(X_right) ... 
            + X_velocity(X_right) - region_x;
        % Do the same thing for the left X boundary (x = 0)
        X_left = X_position + X_velocity < 0;
        X_position(X_left) = X_position(X_left) ... 
            + X_velocity(X_left) + region_x;
        
       % Since we are using logical array, any array position with a '0'
       % means that the position did not go over the canvas boundary.
       % Therefore, we can just simply add the position and velocity
       % for these particles and not worry about the boundary stuff.
       
       % the cool particles are the ones that aren't in x_left or x_right
       the_cool_particles = ~(X_left | X_right);
       X_position(the_cool_particles) = X_position(the_cool_particles) ...
           + X_velocity(the_cool_particles);
       
       % Now we need to worry about the Y boundary...
       % For this, we can check if the Y coordinate is over the limit,
       % and if it is we will just make the velocity negative to turn
       % the particle around!
       % Check if the Y coordinate isn't cool:
       is_Y_not_cool = (Y_position + Y_velocity > region_y ... 
           | Y_position + Y_velocity < 0);
       % For the Y particles that aren't cool (1's in the logical array),
       % we will need to flip their velocities and turn them onto the path
       % of cool.
       Y_velocity(is_Y_not_cool) = -1*Y_velocity(is_Y_not_cool);
       % Now that the particles are all cool, we can update the positions
       Y_position(1,:) = Y_position(1,:) + Y_velocity(1,:);
       
       % Since we're doing a whole bunch of time steps, we're gonna need
       % to save the current timestep results for future plotting.
       X_STORAGE(p,:) = X_position(1,:);
       Y_STORAGE(p,:) = Y_position(1,:);
            
    end
    
    figure(7)
    hist3([X_position',Y_position']);
    view(35,45)
    title('Electron Density Map')
    
    
    figure(8)
    for paths = 1:particles
        subplot(3,1,1);
        plot(X_STORAGE(:,paths),Y_STORAGE(:,paths),'-')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Final Electron Paths')
        hold on
    end
% Box Points for Plotting
    box_pic_x = [0.8e-7 0.8e-7  1.2e-7 1.2e-7];
    box_top_y = [1e-7 0.6e-7 0.6e-7 1e-7];
    box_bottom_y = [0 0.4e-7 0.4e-7 0];
    for row = 1:timestep
        figure(8)
        subplot(3,1,2);
        plot(X_STORAGE(row,:),Y_STORAGE(row,:),'o')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Electrons Flying Around Doing Cooler Things While Having Boxes in the Way')
        hold on
        plot(box_pic_x,box_top_y,'k')
        plot(box_pic_x,box_bottom_y,'k')
        subplot(3,1,3);
        %plot(row,current_temp(row),'.k');
        %title('The Current Temperature Caused by Those Electrons')
        %ylabel('Temperature (K)')
        %xlabel('Time-step')
        %  mean time between collisions is ?mn = 0.2ps
        %legend(['Current Temperature:' num2str(current_temp(row))], ...
        %    ['Avg Time Between Collisions:' num2str(avg_tbc)], ...
        %    ['Mean Free Path:' num2str(MFP)])
        hold on
        xlim([0 timestep])
        ylim([250 350])
        pause(0.00000001)
    end
        
