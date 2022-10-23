% Filename: Problem3
% Simulation of beam deformation with a single external force P=2000N
% with mass spring system 
% The codes are written referencing the lecture notes from Prof.Khalid

%% CITED SOURCES - Cited from Matlab Appendix 
% All the copyrights or the cited sources are reserved to Prof.Khalid
% gradEs.m
% gradEb.m
% hessEs.m
% hessEb.m

%% Given data
clc;
clear;
N=50;      % Number of DOF
dt = 1e-2; % discrete time interval [s]
t = 1;     % simulation time (total time) [s]

% elastic beam
l = 1;        % length [m]
r_o = 0.013;  % outer radius [m]
r_i = 0.011;  % inner radius [m]
E = 70*1e+9;  % elasticity   [Pa]

I = pi/4*(r_o^4-r_i^4);  % moment of Inertia [kg*m^2]
A = pi*(r_o^2-r_i^2);    % area of the sylinder [m^2]

% distance between nodes
dl = l/(N-1); % deltaL

% density [kg/m^3]
rho = 2700; % aluminum

% mass located at each node
m = pi*(r_o^2-r_i^2)*l*rho/(N-1);

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

for i=1:N
    M(2*i-1,2*i-1) = m;          % xi
    M(2*i,2*i) = M(2*i-1,2*i-1); % yi
end

% external force [N]
P = zeros(2*N,1);
numP = round(0.75/dl) + 1; % nth node where external force P is applied on
P(2*numP) = -2000;

% stretching and bending stiffness of a beam
EA = E*A;   % stretching stiffness
EI = E*I;   % bending stiffness

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

% error tolerance value for the loop
e = EI/l^2*1e-3;

% number of times the loop would be repeated
num = round(t/dt);

% boundary condition
fixed = [1;2;2*N]; 
free = 3:2*N-1;

% maximum deformation of the beam along y axis
yMax = zeros(num,1);
yMax(1) = min(q);

for i=2:num
   % initial trial
    q=q0;
    q_free=q(free);

    % initialize error value
    err = 10*e;

    while err>e
        f=M/dt*((q-q0)/dt-u)-P;  % discrete equation of motion
        J=M/dt^2;                % corresponding Jacobian

        for j=1:N-1        
            % elastic force calculation
            % linear spring system for ith sphere and (i+1)th sphere
            xk = q(2*j-1);   % xj
            yk = q(2*j);     % yj
            xkp1 = q(2*j+1); % x_(j+1)
            ykp1 = q(2*j+2); % y_(j+1)
        
            % discrete stretching force
            dF = gradEs(xk, yk, xkp1, ykp1, dl , EA); 
            dJ = hessEs(xk, yk, xkp1, ykp1, dl, EA);
        
            % create index vector
            index = [2*j-1; 2*j; 2*j+1; 2*j+2]; 
        
            % add the calculated values
            f(index) = f(index) + dF;
            J(index,index) = J(index,index) + dJ;
        end

        % bending calculation
        % there is no spring at the first node
        for k=2:N-1
            xkm1 = q(2*k-3); % x_(k-1)
            ykm1 = q(2*k-2); % y_(k-1)
            xk = q(2*k-1);   % xk
            yk = q(2*k);     % yk
            xkp1 = q(2*k+1); % x_(k+1)
            ykp1 = q(2*k+2); % y_(k+1)
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, 0, dl, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, 0, dl, EI); 
            
            % create index vector
            index = [2*k-3; 2*k-2; 2*k-1; 2*k; 2*k+1; 2*k+2];
        
            % add the calculated values
            f(index) = f(index) + dF;
            J(index,index) = J(index,index) + dJ;
        end

        % update the info on free elements
        f_free = f(free);
        J_free = J(free, free);

        % update DOF - Newton Raphson
        q_free = q_free-J_free\f_free;

        % store it back in DOF vector
        q(free) = q_free; 
        
        % update error value
        err = sum(abs(f_free));
    end

    % update velocity
    u = (q-q0)/dt;

    % update old position
    q0 = q;

    % update ymax
    yMax(i) = min(q);
end

% plot the maximum deformation value as a function of time
figure(1);
time = (1:num)*dt;
plot(time, yMax, 'm-');
axis equal
xlabel('time [s]');
ylabel('Ymax [m]');
title('Maximum diformation of the beam as a function of time [m]');





