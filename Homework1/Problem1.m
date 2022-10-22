% Filename: Problem1
% Implicit method simulation of Rigid Spheres and Elastic Beam Falling 
% in Viscous Flow
% The codes are written referencing the lecture notes from Prof.Khalid

%% CITED SOURCES - Cited from Matlab Appendix 
% All the copyrights or the cited sources are reserved to Prof.Khalid
% gradEs.m
% gradEb.m
% hessEs.m
% hessEb.m

%% Implicit method
clc;
clear;

N = 3;        % Number of DOF

g = 9.8;      % gravitational acceleration [m/s^2]
E = 1e+9;     % Young's modulus

% radius of spheres [m]
r1 = 0.005; 
r2 = 0.025; 
r3 = 0.005;

% density [kg/m^3]
rho_m = 7000; % sphere
rho_f = 1000; % fluid

% fluid viscosity [Pa-s]
mu = 1000;

% elastic beam
l = 0.1;      % length [m]
r0 = 0.001;   % cross sectional radius [m]

% distance between nodes
dl = l/(N-1); % deltaL

% stretching and bending stiffness of a beam
EA = E*pi*r0^2;     % stretching stiffness
EI = E*pi*r0^4/4;   % bending stiffness

% total time [s]
t = 10; 

% discrete time interval [s]
dt = 1e-2; % implicit method

% initial position
% update the initial positions - x and y positions - of each sphere
positions = zeros(N,2);

for i=1:N
    positions(i,1)=(i-1)*dl;
    positions(i,2)=0;
end

% mass
% update the mass of each sphere
M = zeros(2*N,2*N);

% mass of each sphere
m1 = 4*pi/3*r1^3*rho_m;
m2 = 4*pi/3*r2^3*rho_m;
m3 = 4*pi/3*r3^3*rho_m;

M(1,1)=m1; % x1
M(2,2)=m1; % y1
M(3,3)=m2; % x2
M(4,4)=m2; % y2
M(5,5)=m3; % x3
M(6,6)=m3; % y3

% viscous damping matrix
% update viscous damping of each sphere
C = zeros(2*N,2*N);

% viscous damping of each sphere
c1 = 6*pi*mu*r1;
c2 = 6*pi*mu*r2;
c3 = 6*pi*mu*r3;

C(1,1) = c1;
C(2,2) = c1;
C(3,3) = c2;
C(4,4) = c2;
C(5,5) = c3;
C(6,6) = c3;

% weight
% update weight vector elements
W = zeros(2*N,1);

rho = rho_m-rho_f; % density

% weight of each sphere
w1 = -4/3*pi*r1^3*rho*g;
w2 = -4/3*pi*r2^3*rho*g;
w3 = -4/3*pi*r3^3*rho*g;

W(2,1) = w1;
W(4,1) = w2;
W(6,1) = w3;

% initial DOF vector
q0 = zeros(2*N,1);

% imply the position information in position matrix to this new matrix
for i=1:N
    q0(2*i-1)= positions(i,1); % x position
    q0(2*i) = positions(i,2);  % y position
end

% initialize new DOF and velocity vector with initial conditions
q=q0;
u=(q-q0)/dt;

% tolerance error value for the loop for Newton Raphson method
e = EI/l^2*1e-3;

% number of times the loop would be repeated
num = round(t/dt);

% middle node position and velocity along y axis
x2 = zeros(num,1);
v2 = zeros(num,1);
midN = (N+1)/2;    % middle node
x2(1) = u(2*midN); % update initial position
v2(1) = u(2*midN); % update initial velocity

for i=2:num
    % initial trial
    q=q0;

    % initialize error value
    err = 10*e;

    while err>e
        f=M/dt*((q-q0)/dt-u);    % discrete equation of motion
        J=M/dt^2;                % corresponding Jacobian

        % elastic force calculation
        % linear spring system for sphere 1 and sphere 2
        xk = q(1);   % x1
        yk = q(2);   % y1
        xkp1 = q(3); % x2
        ykp1 = q(4); % y2

        % discrete stretching force
        dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, dl, EA); % corresponding Jacobian
        
        % add the calculated values
        f(1:4) = f(1:4)+ dF;
        J(1:4, 1:4) = J(1:4, 1:4) + dJ;

        % linear spring system for sphere 2 and sphere 3
        xk = q(3);   % x3
        yk = q(4);   % y3
        xkp1 = q(5); % x4
        ykp1 = q(6); % y4

        % discrete stretching force
        dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, dl, EA); % corresponding Jacobian

        % add the calculated values
        f(3:6) = f(3:6)+ dF;
        J(3:6, 3:6) = J(3:6, 3:6) + dJ;

        % bending force calculation
        xkm1 = q(1); % x1
        ykm1 = q(2); % y1
        xk = q(3);   % x2
        yk = q(4);   % y2
        xkp1 = q(5); % x3
        ykp1 = q(6); % y3
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, 0, dl, EI);
        dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, 0, dl, EI);

        % add the calculated values
        f(1:6) = f(1:6)+ dF;
        J(1:6, 1:6) = J(1:6, 1:6) + dJ;

        % viscous force calculation
        % discrete viscous force
        f = f+C*(q-q0)/dt;
        J=J+C/dt; % corresponding Jacobian

        % weight calculation
        f = f-W;

        % update DOF - Newton raphson
        q = q-J\f;
        
        % update error value
        err = sum(abs(f));
    end

    % update velocity
    u = (q-q0)/dt;

    % update old position
    q0 = q;

    % update middle node position
    x2(i) = q(2*midN);

    % update middle node velocity
    v2(i) = u(2*midN);
end

% plot the shape of structure
xpos = q(1:2:end);
ypos = q(2:2:end);
figure(1);
plot(xpos, ypos, 'mo-');
axis equal
xlabel('x [m]');
ylabel('y [m]');
title('Shape of the structure');

% plot the position of the sphere in the middle node
figure(2);
time = (1:num)*dt;
plot(time, x2, 'm-');
xlabel('Time [s]');
ylabel('Position of the sphere in middle node [m]');
title('Position of the sphere in middle node along y axis');

% plot the velocity of the sphere in the middle node
figure(3);
time = (1:num)*dt;
plot(time, v2, 'm-');
xlabel('Time [s]');
ylabel('Velocity of  the sphere in middle node [m/s]');
title('Velocity of the sphere in middle node along y axis');

%% Explicit method
clc;
clear;

N = 3;        % Number of DOF

g = 9.8;      % gravitational acceleration [m/s^2]
E = 1e+9;     % Young's modulus

% radius of spheres [m]
r1 = 0.005; 
r2 = 0.025; 
r3 = 0.005;

% density [kg/m^3]
rho_m = 7000; % sphere
rho_f = 1000; % fluid

% fluid viscosity [Pa-s]
mu = 1000;

% elastic beam
l = 0.1;      % length [m]
r0 = 0.001;   % cross sectional radius [m]

% distance between nodes
dl = l/(N-1); % deltaL

% stretching and bending stiffness of a beam
EA = E*pi*r0^2;     % stretching stiffness
EI = E*pi*r0^4/4;   % bending stiffness

% simulation time (total time) [s]
t = 10; 

% discrete time interval [s]
dt = 1e-5; % explicit method

% initial position
% update the initial positions - x and y positions - of each sphere
positions = zeros(N,2);

for i=1:N
    positions(i,1)=(i-1)*dl;
    positions(i,2)=0;
end

% mass
% update the mass of each node
M = zeros(2*N,2*N);

% mass of each sphere
m1 = 4*pi/3*r1^3*rho_m;
m2 = 4*pi/3*r2^3*rho_m;
m3 = 4*pi/3*r3^3*rho_m;

M(1,1)=m1; % x1
M(2,2)=m1; % y1
M(3,3)=m2; % x2
M(4,4)=m2; % y2
M(5,5)=m3; % x3
M(6,6)=m3; % y3

% viscous damping matrix
% update viscous damping of each sphere
C = zeros(2*N,2*N);

% viscous damping of each sphere
c1 = 6*pi*mu*r1;
c2 = 6*pi*mu*r2;
c3 = 6*pi*mu*r3;

C(1,1) = c1;
C(2,2) = c1;
C(3,3) = c2;
C(4,4) = c2;
C(5,5) = c3;
C(6,6) = c3;

% weight
% update weight vector elements
W = zeros(2*N,1);

rho = rho_m-rho_f; % density

% weight of each sphere
w1 = -4/3*pi*r1^3*rho*g;
w2 = -4/3*pi*r2^3*rho*g;
w3 = -4/3*pi*r3^3*rho*g;

W(2,1) = w1;
W(4,1) = w2;
W(6,1) = w3;

% initial DOF vector
q0 = zeros(2*N,1);

% imply the position information in position matrix to this new matrix
for i=1:N
    q0(2*i-1)= positions(i,1); % x position
    q0(2*i) = positions(i,2);  % y position
end

% initialize new DOF and velocity vector with initial conditions
% initialize force matrix as well
q=q0;
u = (q-q0)/dt;
f = zeros(2*N,1);

% number of times the loop would be repeated
num = round(t/dt);

% middle node position and velocity along y axis
x2 = zeros(num,1);
v2 = zeros(num,1);
midN = (N+1)/2;    % middle node
x2(1) = u(2*midN); % update initial position
v2(1) = u(2*midN); % update initial velocity

for i=2:num   
    % update force
    f = zeros(2*N,1);

    % elastic force calculation
    % linear spring system for sphere 1 and sphere 2
    xk = q0(1);   % x1
    yk = q0(2);   % y1
    xkp1 = q0(3); % x2
    ykp1 = q0(4); % y2

    % discrete stretching force
    dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
        
    % add the calculated values
    f(1:4) = f(1:4)+ dF;

    % linear spring system for sphere 2 and sphere 3
    xk = q0(3);   % x3
    yk = q0(4);   % y3
    xkp1 = q0(5); % x4
    ykp1 = q0(6); % y4

    % discrete stretching force
    dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);

    % add the calculated values
    f(3:6) = f(3:6)+ dF;

    % bending force calculation
    xkm1 = q0(1); % x1
    ykm1 = q0(2); % y1
    xk = q0(3);   % x2
    yk = q0(4);   % y2
    xkp1 = q0(5); % x3
    ykp1 = q0(6); % y3
    dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, 0, dl, EI);

    % add the calculated values
    f(1:6) = f(1:6)+ dF;

    % viscous force calculation
    % discrete viscous force
    f = f+C*u;

    % weight calculation
    f = f-W;
    
    % update new position
    for j=1:6
        q(j) = q0(j)+dt*(u(j)-dt/M(j,j)*f(j));
    end

    % update velocity
    u = (q-q0)/dt;

    % update old position
    q0 = q;

    % update middle node position
    x2(i) = q(2*midN);

    % update middle node velocity
    v2(i) = u(2*midN);
end

% plot the shape of structure
xpos = q(1:2:end);
ypos = q(2:2:end);
figure(4);
plot(xpos, ypos, 'mo-');
axis equal
xlabel('x [m]');
ylabel('y [m]');
title('Shape of the structure');





