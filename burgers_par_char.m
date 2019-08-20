%% ****** Domain ******
xStart = 0; xEnd = 1; % Range of Domain
nb_cells = 100;  % Number of cells
nb_particles = 100; % Number of particles per cell
x = linspace(xStart, xEnd, nb_cells + 1); % discretize into cells
dx = x(2) - x(1); % grid size
%% ****** Particles ******
total_nb_particles = nb_particles * nb_cells;
par_old = zeros(3, total_nb_particles); % particle history -old
par_new = par_old; % particle history - new
%% ****** Cell centre calculation ******
cell_centre = zeros(1, nb_cells);
for i = 1:nb_cells
    cell_centre(i) = x(i) + (dx / 2);
end
%% ****** Memory allocation for cell velocity *******
cell_vel = 0 .* cell_centre;  % old cell velocity
new_cell_vel = 0 .* cell_vel; % New cell veloctiy
particles = 0 .* cell_centre; % current cell numbers of particles
%% ****** Cell velocity initialization ********
cell_vel=sin(2*pi*cell_centre);
%%
CFL = 0.1;
final_time = 0.14; % final time
dt = CFL * dx / max(cell_vel)
nt = uint32(final_time/dt)
%%  ******* Particle position initialization *******
rng('default');
rng(1);
for i = 1:nb_cells
    for j = (((i - 1) * nb_particles) + 1):(i * nb_particles)
        par_old(1, j) = x(i) + dx * rand(1, 1); % particle position
    end
end
%% ****** Partical velocity initialization *******
for i = 1:nb_cells
    for j = (((i - 1) * nb_particles) + 1):(i * nb_particles)
        par_old(2, j) = cell_vel(1, i); % particle velocity
        par_old(3, j) = i; % particle cell number
    end
end
% Copy the current particle history to the another avarible for maipulation
par_new = par_old;
%% main loop
fileID1=fopen('U.txt','w');
for time = 1:nt
    for i = 1:total_nb_particles
        par_new(1, i) = par_old(1, i) + (dt * par_old(2,i));
        if (par_new(1, i) > xEnd)
            par_new(1, i) = par_new(1, i) - xEnd;
        elseif (par_new(1, i) < xStart)
            par_new(1, i) = par_new(1, i) + xEnd;
        end
    end
    for i = 1:total_nb_particles
        position_x = par_new(1, i);
        xl = 1;
        xr = nb_cells + 1;
        xm = floor((xl + xr) / 2);
        [xl, xr] = find_cell_x(position_x, xl, xr, xm, x);
        par_new(3, i) = xl;
    end
    new_cell_vel = zeros(1, nb_cells);
    particles = zeros(1, nb_cells);
    for i = 1:total_nb_particles
        c = par_new(3, i);
        new_cell_vel(1, c) = new_cell_vel(1, c) + par_new(2, i);
        particles(1, c) = particles(1, c) + 1;
    end
    for i = 1:nb_cells
        if (particles(1, i) == 0)
            new_cell_vel(1, i) = 0;
        else
            new_cell_vel(1, i) = new_cell_vel(1, i) / particles(1, i);
        end
    end
 
    for i = 1:total_nb_particles
        c = par_new(3, i);
        par_new(2, i) = new_cell_vel(1, c);
    end
    par_old = par_new;
    fprintf(fileID1,'%16f, ',new_cell_vel);
    fprintf(fileID1,'\n');
end
%% EACT SOLUTION using characteristics
X0=cell_centre;  % initial position of characteristics
u = cell_vel; % initial velocity 
counter =0;
fileID2=fopen('u.txt','w');
for j = 1:nt  % time loop
    counter = counter+1;
    X = X0 - (u .* (counter*dt) ); % position update -> X(0)= X(t)-u*dt
    for i=1:nb_cells
        if (X(1, i) > xEnd)
            X(1, i) = x(1, i) - xEnd;
        elseif (X(1, i) < xStart)
            X(1, i) = X(1, i) + xEnd;
        end
    end
    u=sin(2*pi*X);
    fprintf(fileID2,'%16f, ',u);
    fprintf(fileID2,'\n');
end
error = norm(new_cell_vel-u)
fclose(fileID1);
fclose(fileID2);
fileID3=fopen('e.txt','w');
fprintf(fileID3,'%16f\n',error);
fclose(fileID3);
