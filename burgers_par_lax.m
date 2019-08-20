% GAUSSIAN INITIAL CONDITION
%u0 = exp(-200*(x-0.25).^2);
%u_exact = exp(-200*(x-0.25-a*(dt*j)).^2);
% WORKS ALSO FOR SIN INITIAL CONDITION

%% this is inviscid burgers case- DON'T CHANGE N=100, par=50
final_time = 0.1; t=0;
xStart = 0; xEnd = 1; nb_cells = 1000; nb_particles = 100;
x = linspace(xStart, xEnd, nb_cells + 1);
dx = x(2) - x(1);
cell_centre = zeros(1, nb_cells);
for i = 1:nb_cells
    cell_centre(i) = x(i) + (dx / 2);
end

cell_vel = 0 .* cell_centre;
new_cell_vel = 0 .* cell_vel;

cell_vel = sin(2 * pi * cell_centre); %CHANGE VELOCITY HERE
CFL=1;
dt = CFL*dx/max(cell_vel);
nt = uint32(final_time / dt);
total_nb_particles = nb_particles * nb_cells;
par_old = zeros(4, total_nb_particles);
par_new = par_old;
U = zeros(nt + 1, nb_cells);
%% Initial uniform distribution of articles
rng('default');
rng(1);
for i = 1:nb_cells
    for j = (((i - 1) * nb_particles) + 1):(i * nb_particles)
        par_old(1, j) = x(i) + dx * rand(1, 1);
    end
end
%% cell velocity initialization
U(1, :) = cell_vel(1, :);
%% partical velocity initialization
for i = 1:nb_cells
    for j = (((i - 1) * nb_particles) + 1):(i * nb_particles)
        par_old(2, j) = cell_vel(1, i);
        par_old(3, j) = i;
    end
end
par_new = par_old;
%% main loop
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
    U(time + 1, :) = new_cell_vel(1, :);
end
%%
for j=1:nt
    plot(cell_centre,U(j,:),'k.');
    pause(dt)
end
%hold on 
%******************************************
hlf = plot(cell_centre, g(cell_centre), 'b.');
ue =zeros(nt+1,nb_cells);
u = cell_vel;
utemp = u;
ue(1,:)=u;
for time=1:nt
%while (t+dt<final_time)
    % Lax-Friedrichs
    for i=2:nb_cells-1
        dflux = 0.5*u(i+1)^2 - 0.5*u(i-1)^2;
        utemp(i) = 0.5*(u(i+1) + u(i-1)) - 0.5*dt/dx* dflux;
    end
    utemp(1) = utemp(nb_cells-1);
    utemp(nb_cells) = utemp(2);

    u = utemp;
    ue(time+1,:)=u;
    %t = t + dt;
    %dt = CFL*dx/max(u);
    set(hlf,'YData',u);
    %plot(cell_centre,u,'b.')
    drawnow;
end
%legend('Char.','LF','Location','northwest');
%%
figure(2)
errors=zeros(1,nt+1);
counter =1:nt+1;
for i=1:nt+1
errors(i) = norm(U(i,:)-ue(i,:));
plot(counter,errors);
%counter =counter+1;
end

%%
function y =g(x)
 %y=x;
y =sin(2*pi*x);
end
