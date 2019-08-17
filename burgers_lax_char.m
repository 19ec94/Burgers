nb_cells = 10000; % number of nodes in x
CFL = 0.99; % Courant number
final_time = 0.1; % final time
% initialization
t = 0;
xStart=0; xEnd=1;
x = linspace(xStart,xEnd,nb_cells+1);
dx = x(2)-x(1);
cell_centre = x + (dx/2);
cell_centre(nb_cells+1)=[];
u= sin(2*pi*cell_centre); % CHANGE VELOCITY HERE FOR dt
dt = CFL*dx/max(u);
nt = uint32(final_time/dt);
%******************************************
%hlf = plot(cell_centre, g(cell_centre), 'b.');
ue =zeros(nt+1,nb_cells);
utemp = u;
ue(1,:)=u;
for time=1:nt
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
    %set(hlf,'YData',u);
    %plot(cell_centre,u,'b.')
    %drawnow;
end
%% Exact char
u_exact = zeros(nt + 1, nb_cells); % stores  exact data for all time   
%u_exact(1,:) = cell_centre; % Initial velocity
u_exact(1,:) = sin(2*pi*cell_centre); % Initial velocity
X0=cell_centre;  % initial position of characteristics
u_E = u_exact(1,:); % initial velocity 
counter =0;
for j = 1:nt  % time loop
    counter = counter+1;
    X = X0 - (u_E .* (counter*dt) ); % position update -> X(0)= X(t)-u*dt
    if (X(1, i) > xEnd)
        X(1, i) = x(1, i) - xEnd;
    elseif (X(1, i) < xStart)
        X(1, i) = X(1, i) + xEnd;
    end
    u_E = sin(2*pi*X);  % solution at the new time  u(X,t) = u((X(0)),t)
    u_exact(j + 1, :) = u_E;% for data analysis
end
%% Error between particle solution and characteristic solution
figure(2)
errors=zeros(1,nt+1);
counter =1:nt+1;
for i=1:nt+1
errors(i) = norm(u_exact(i,:)-ue(i,:));
plot(counter,errors);
%counter =counter+1;
end
%Error = norm(u_exact(nt+1,:)-ue(nt+1,:));
%%
% plot(cell_centre,u_exact(nt+1,:),'*',cell_centre,u,'o');
function y =g(x)
y =sin(2*pi*x);
end
