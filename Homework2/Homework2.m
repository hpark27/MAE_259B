% Filename: Homework2
% 3D simulation of Discrete Elastic Rod under gravity for 5 seconds
% The codes are written referencing the lecture notes from Prof.Khalid

%% CITED SOURCES - Cited from Matlab Appendix 
% All the copyrights or the cited sources are reserved to Prof.Khalid
% computekappa.m
% crossMat.m
% gradEb_hessEb.m
% gradEs_hessEs.m
% gradEt_hessEt.m
% signedAngle.m
% plotRod.m

%%
clear;
clc;

% global variable
global Fg;       % gravitational force
global mMat;     % mass matrix
global dt;       % time stamp
global kappaBar;

% for bending
global EI;       % bending stiffness
global vorL;     % voroni length

% for twist
global GJ;       % twisting stiffness

% for stretch
global EA;       % stretching stiffness
global refL;     % reference length

% input
node = 50;       % number of nodes
edge = node-1;   % number of edges
N = 3*node+edge; % number of DOF

dt = 0.01;       % time stamp

% rod parameters
l = 0.2;         % rod length
nR = 0.02;       % natural radius
r0 = 0.001;      % cross sectional radius

rho = 1000;      % density [kg/m^3]

Y = 10*1e+6;     % Young's modulus
v = 0.5;         % Poisson's ratio
G = Y/(2*(1+v)); % Shear modulus

g = [0;0;-9.8];  % gravitational acceleration [m/s^2] only in z direction

t = 5;           % total simulation time

% stiffness variables
EI = Y*pi*r0^4/4; % Bending
GJ = G*pi*r0^4/2; % Shearing
EA = Y*pi*r0^2;   % Stretching

% Error tolerance
tol = EI/l^2*1e-6;

m = pi*r0^2*l*rho; % total mass of rod
dm = m/edge;       % mass for each edge

% mass matrix
M = zeros(N,1);

% update mass matrix elements by looping over nodes
for i = 1:node
    index = [4*i-3; 4*i-2; 4*i-1]; 

    if i == 1
        M(index) = dm/2; 
    elseif i == node
        M(index) = dm/2;
    else
        M(index) = dm;
    end
end

% update mass matrix elements by looping over edges
for i = 1:edge
    M(4*i) = 1/2*dm*r0^2;
end

mMat = diag(M); % set up diagonal matrix using mass matrix

nodes = zeros(node, 3); % initialize three dimensional node matrix
dTheta = l/(nR*edge);   % delta theta

% update node matrix elements
for i = 1:node
    nodes(i,1) = nR*cos((i-1)*dTheta); % x coordinate of c_th node
    nodes(i,2) = nR*sin((i-1)*dTheta); % y coordinate of c_th node
    nodes(i,3) = 0;                    % z coordinate of c_th node
end
 
% create initial DOF vector
q0 = zeros(N, 1); % ndof = 4N-1

% update DOF vector elements by looping over nodes
for i = 1:node
    index = [4*i-3; 4*i-2; 4*i-1];
    q0(index) = nodes(i,:); % update array element
end

% create reference length - the length of edge - vector
refL = zeros(edge,1);

% update vector elements by looping over edges
for i = 1:edge
    dx = nodes(i+1,:)-nodes(i,:);
    refL(i) = norm(dx);
end

% create a vector to for Voroni length - length associated with each node
% loop over nodes
vorL = zeros(node,1);

% update vector elements by looping over nodes
for i = 1:node
    if i == 1
        vorL(i) = 0.5*refL(i);
    elseif i == node
        vorL(i) = 0.5*refL(i-1);
    else
        vorL(i) = 0.5*(refL(i-1)+refL(i));
    end
end

% create a vector for gravity
Fg = zeros(N,1);

% update vector elements by looping over node
for i=1:node
    index = [4*i-3; 4*i-2; 4*i-1];
    
    % gravitational force = m*g
    Fg(index) = M(index).*g;
end

% reference frames
% space parallel transport at t=0
ref1 = zeros(edge, 3);        % 1st reference director
ref2 = zeros(edge, 3);        % 2nd reference director
tangent = computeTangent(q0); % Tangent

tan0 = tangent(1,:);          % tangent on the first edge
t0 = [0;0;-1];                % arbitrary vector
a1Tmp = cross(tan0, t0);      % cross product of the two vector

if abs(a1Tmp) < 1e-6 % when a1Tmp is really small, and it is close to 0
    t0 = [0;1;0];
    a1Tmp = cross(tan0, t0);
end

% plug a1tmp into frame after normalize it
ref1(1,:) = a1Tmp /norm(a1Tmp);
ref2(1,:) = cross(tangent(1,:), ref1(1,:));

% space parallel transport to construct the reference frame
% since the computations for the first edge is done, start from the
% second edge
for i = 2:edge
    itan0 = tangent(i-1,:);   % tangent on (i-1)th edge
    itan1 = tangent(i,:);     % tangent on ith edge

    ref1_tan0 = ref1(i-1, :); % first reference frame from c-1 edge
    % transport it from previous edge to next edge
    ref1_tan1 = parallel_transport(ref1_tan0, itan0, itan1); 
    
    % store the value after normalizing it
    ref1(i,:) = ref1_tan1/norm(ref1_tan1); 
    ref2(i,:) = cross(itan1,ref1(i,:));
end

% material frame
theta = q0(4:4:end); % angle vector
% rotate reference frame ref1 and ref2 by angle theta
[m1, m2] = computeMaterialDirectors(ref1, ref2, theta);

% reference twist
% it is zero beccause reference frames are initialized using space parallel
% transport
refTwist = zeros(node,1);

% natural curvature
kappaBar = getkappa(q0, m1, m2);

% fixed and free DOF
fixed = 1:7; % fixed DOF
free = 8:N;  % free DOf

% time stepping scheme
nStep = round(t/dt);   % number of steps
cTime = 0;             % current time
endZ = zeros(nStep,1); % z-coordinate of last node

% initialize position and velocity
q = q0;             % position  
u = zeros(size(q)); % velocity

% timestep
for i=1:nStep
    [q, u, ref1, ref2] = objfun(q0, u, ref1, ref2, free, tol, refTwist);

    % update current time
    cTime = cTime + dt;
    
    % update position
    q0 = q;

    % store data
    endZ(i) = q(end);

    % plot
    if mod(i, 100) == 0
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(ref1, ref2, theta); 
        plotrod(q, ref1, ref2, m1, m2, cTime);
    end
end

% plot
figure(2);
timearray = (1:1:nStep)*dt;
plot(timearray,endZ,'mo-');
box on
xlabel('Time [s]');
ylabel('Z coordinate of last node [m]');




















