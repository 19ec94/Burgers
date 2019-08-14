%% ****** Time ******
dt = 0.001; % time step size
final_time = 0.3; % final time
nt = uint32(final_time / dt); % total number of time steps
%% ****** Domain ******
xStart = 0; xEnd = 1; % Range of Domain
nb_cells = 100;  % Number of cells
nb_particles = 30; % Number of particles per cell
x = linspace(xStart, xEnd, nb_cells + 1); % discretize into cells
dx = x(2) - x(1); % grid size
%% ****** Particles ******
total_nb_particles = nb_particles * nb_cells;
par_old = zeros(3, total_nb_particles); % particle history -old
par_new = par_old; % particle history - new
% par_old(1,:) -> stores particles positions 
% par_old(2,:) -> stores particles velocities
% par_old(3,:) -> Current cell numbers of the particles
%% ****** Cell centre calculation ******
cell_centre = zeros(1, nb_cells);
for i = 1:nb_cells
    cell_centre(i) = x(i) + (dx / 2);
end
%% ****** Memory allocation for cell velocity *******
cell_vel = 0 .* cell_centre;  % old cell velocity
new_cell_vel = 0 .* cell_vel; % New cell veloctiy
particles = 0 .* cell_centre; % current cell numbers of particles
% To store all velocity results for all time for data analysis
U = zeros(nt + 1, nb_cells);
%%  ******* Particle position initialization *******
% Particles in each cell are distributed uniformly. This is why we use rand
% function here. we use the rng function for debugging purpose. This
% function ensures we get the same random number for every repitition of
% simuation
rng('default');
rng(1);
for i = 1:nb_cells
    for j = (((i - 1) * nb_particles) + 1):(i * nb_particles)
        par_old(1, j) = x(i) + dx * rand(1, 1); % particle position
    end
end

%% ****** Cell velocity initialization ********
% We evaluate cell velocity at the centre point of each cell. And this will
% be our cell velocity
cell_vel =cell_centre;  % Linear Initial condition
U(1, :) = cell_vel(1, :); % store it in our variable for data analysis

%% ****** Partical velocity initialization *******
% We assign each particle a velocity. Each particle will get a velocity of
% the cells in which it is currently.  For example 5th particle in 6th cell
% get a velocity of 6th cell. par old(1,5) = cell vel(1,6)
for i = 1:nb_cells
    for j = (((i - 1) * nb_particles) + 1):(i * nb_particles)
        par_old(2, j) = cell_vel(1, i); % particle velocity
        par_old(3, j) = i; % particle cell number
    end
end
% Copy the current particle history to the another avarible for maipulation
par_new = par_old;
%% main loop
for time = 1:nt
    % Position update
    for i = 1:total_nb_particles
        par_new(1, i) = par_old(1, i) + (dt * par_old(2,i)); %xNew=xOld+dt*a
    end
    %Find where are my particles are and store it in par_old(3,:)    
    for i=1:total_nb_particles
        for j=1:nb_cells  
            if( par_new(1,i) >= x(1,j)) && ( par_new(1,i) <= x(1,j+1))
                par_new(3,i)=j;
            end
        end
    end
    for i=1:total_nb_particles
            if( par_new(1,i) > xEnd) || ( par_new(1,i) < xStart)
                par_new(3,i)=-1;
            end
    end   
    % Calculate total velocity for each cell
    new_cell_vel = 0 .* new_cell_vel;
    particles = 0 .* particles;
    for i = 1:total_nb_particles
        c = par_new(3, i);
        if (c == -1)
            % do nothing
        else
            new_cell_vel(1, c) = new_cell_vel(1, c) + par_new(2, i);
            particles(1, c) = particles(1, c) + 1;
        end
    end
    % divide each cell by number of current particles
    for i = 1:nb_cells
        if (particles(1, i) == 0)
            new_cell_vel(1, i) = 0;
        else
            new_cell_vel(1, i) = new_cell_vel(1, i) / particles(1, i);
        end
    end
    % distribute new velocity to each particles
    for i = 1:total_nb_particles
        c = par_new(3, i);
        if (c==-1)
            par_new(2,i) = 0;% Mute the particles that are out of boundary
        else
            par_new(2, i) = new_cell_vel(1, c);
        end
    end
    par_old = par_new; % copy old data
    U(time + 1, :) = new_cell_vel(1, :);% This is for data analysis
end
%% EACT SOLUTION using characteristics
u_exact = zeros(nt + 1, nb_cells); % stores  exact data for all time   
% cell_centre is our X0;
u_exact(1,:) = cell_centre; % Initial velocity
X=cell_centre;  % initial position of characteristics
u = u_exact(1,:); % variable for manipuation 
for j = 1:nt
    X = X - (u .* (dt) ); % position update -> X(t)= X(t-1)-u*dt
    u = X;  % solution at the new time  u(X,t) = u((X0-u*t),t)
    u_exact(j + 1, :) = u; % for data analysis
end
%% Error between particle solution and characteristic solution
Error = norm(U(nt+1,:)-u_exact(nt+1,:));
%%
%% Plot 
% If you want to plot, delete "%{" in the beginning and "%}" end of the
% following code snippet

counter=0;
for j=1:nt
    counter =counter+1;
    t = counter*dt;
    plot(cell_centre,cell_centre,cell_centre,U(j,:),cell_centre,u_exact(j,:));
    legend('initial','particle','exact');
    title(['t =',num2str(t)]);
    pause(dt)
end
