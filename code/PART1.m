   % PART 1: Monte Carlo Simulator from Assignment 1 with Adjustments

    particles = 50;
    
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

    
    %% PAINTING THE CANVAS ------------------------------------------------
    
    % Go through loop for specified number of time steps
    for p = 1:timestep
        
        % Apply constant acceleration due to E field
        X_velocity(:,:) = X_velocity(:,:) + vel_fromA;
        
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
       

       % Calculate current for the plots
       vel = sqrt((X_velocity/v_change).^2 +(Y_velocity/v_change).^2);
       avg_velocity = sum(vel)/particles;
       mu = (avg_velocity)/E_FIELD;
       CURRENT_STORAGE(p) = C.q_0*elecC*mu*E_FIELD/(region_x*region_y);
       
       % Finally, calculate temperature
       current_temp = 0.5*C.m_n*vel.^2 / (2*C.kb);
       average_temp = sum(current_temp)/particles;
       % Mean free path calculation
       MFP = sum(mfp_vec)/particles;
       avg_tbc = MFP/avg_velocity;
       
    end
    
    %% SHOWING THE ART ----------------------------------------------------
    % Everything is calculated! Now let's start plotting...
    disp('The electric field is calculated to be:')
    disp(E_FIELD)
    disp('The force due to this field is:')
    disp(E_FORCE)
    disp('The acceleration due to this force is:')
    disp(acceleration)
    
    % Current Plot
    disp('The current plot below is calculated using current = q_0 * electron concentration * mu * electric field / area')
    figure(1)
    plot(linspace(1,timestep,timestep),CURRENT_STORAGE)
    title('Current Plot over Time')
	xlabel('Time Step')
    ylabel('Current (A/cm^2)')
    
    % --------------------------------------------------------------------
    
    % Temperature Map
    % Set up the X and Y linspace for the meshgrid
    X_canvas = linspace(0,region_x,20);
    Y_canvas = linspace(0,region_y,20);

    % Group x and y position values into bins defined in x and y canvas
    x_bin = discretize(X_position,X_canvas);
    y_bin = discretize(Y_position,Y_canvas);

    % Set up the temperature bins 
    TEMP_BINS = zeros(20,20);
    % Start going through the x and y values
    for x = 1:20
        for y=1:20
            % Define the x and y bins as logic values to easily assign
            % every value in the x and y velocities to a proper bin
            logicX = x_bin == x;
            logicY = y_bin == y;
            logic = logicX & logicY;
            % Whichever values in X and Y velocity fit into the correct bin
            % Will be added to be summed
            sum_x = sum(X_velocity(logic)) / v_change;
            sum_y = sum(Y_velocity(logic)) / v_change;
            
            velocity_avg = sqrt((sum_x)^2 + (sum_y)^2);
            % Use average velocity to calculate temperature
            % and put it into the proper bin
            TEMP_BINS(x,y) = 0.5*C.m_n*velocity_avg.^2 / (2*C.kb);
        end
    end
    
    % Plot the bins as a temperature map
    figure(2)
    surf(TEMP_BINS)
    title('Electron Temperature Map')
    colorbar
    view(2)
    
    % --------------------------------------------------------------------
    
    % Current Density Map
    % Set up the X and Y linspace for the meshgrid
    X_canvas = linspace(0,region_x,20);
    Y_canvas = linspace(0,region_y,20);

    % Group x and y position values into bins defined in x and y canvas
    x_bin = discretize(X_position,X_canvas);
    y_bin = discretize(Y_position,Y_canvas);

    % Set up the temperature bins 
    DENSITY_BINS = zeros(20,20);
    
    for x = 1:20
        for y=1:20
            % Define the x and y bins as logic values to easily assign
            % every value in the x and y velocities to a proper bin
            logicX = x_bin == x;
            logicY = y_bin == y;
            logic = logicX & logicY;
            DENSITY_BINS(x,y) = sum(x_bin(logic))/x + sum(y_bin(logic))/y;
        end
    end

    figure(3)
    surf(DENSITY_BINS)
    title('Current Density Map')
    colorbar
    view(2)
            
    
    
    figure(4)
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
    for row = 1:timestep
        figure(4)
        subplot(3,1,2);
        plot(X_STORAGE(row,:),Y_STORAGE(row,:),'o')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Electrons Flying Around Doing Cool Things inside of an E Field')
        hold on
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
        pause(0.0001)
    end
    

    
    