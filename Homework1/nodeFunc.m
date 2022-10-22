% a function to run an implicit method simulation with difference node
% number input
% input: number of nodes output: terminal velocity of the system
function n1 = nodeFunc(node)
    N = node;  % Number of DOF
    dt = 1e-2; % discrete time interval [s]
    t = 50;    % simulation time (total time) [s]
    
    % elastic beam
    l = 0.1;      % length [m]
    r0 = 0.001;   % cross sectional radius [m]
    
    % distance between nodes
    dl = l/(N-1); % deltaL
    
    % radius of spheres
    r_mid = 0.025; % middle sphere
    r_s = dl/10;   % remaining spheres on nodes
    
    % middle node index calculation
    midN = (N+1)/2;
    
    g = 9.8;      % gravitational acceleration [m/s^2]
    E = 1e+9;     % Young's modulus
    
    % density [kg/m^3]
    rho_m = 7000; % sphere
    rho_f = 1000; % fluid
    
    % fluid viscosity [Pa-s]
    mu = 1000;
    
    % stretching and bending stiffness of a beam
    EA = E*pi*r0^2;     % stretching stiffness
    EI = E*pi*r0^4/4;   % bending stiffness
    
    % set up radius vector
    r = zeros(N,1);
    r(:)=r_s;        % radius of spheres in nodes
    r(midN) = r_mid; % radius of sphere in middle
   
    % initial position
    % update the initial positions - x and y positions - of each sphere
    positions = zeros(N,2);
    
    for i=1:N
        positions(i,1)=(i-1)*dl;
        positions(i,2)=0;
    end
    
    % mass
    % update mass of each node
    M = zeros(2*N,2*N);
    
    for i=1:N
        % mass of a sphere at ith node
        mi = 4*pi/3*r(i)^3*rho_m;
    
        M(2*i-1,2*i-1) = mi;         % xi
        M(2*i,2*i) = M(2*i-1,2*i-1); % yi
    end
    
    % viscous damping
    % update viscous damping on each node
    C = zeros(2*N,2*N);
    
    for i=1:N
        % viscous damping on a sphere at ith node
        ci = 6*pi*mu*r(i);
    
        C(2*i-1,2*i-1) = ci;         % xi
        C(2*i,2*i) = C(2*i-1,2*i-1); % yi
    end
    
    % weight
    % update weight of each node
    W = zeros(2*N,1); 
    
    rho = rho_m - rho_f; % density
    
    for i=1:N
        % weight of a sphere at ith node
        wi = -4*pi*r(i)^3/3*rho*g;
        W(2*i-1) =  0; % xi
        W(2*i) = wi;   % yi 
    end
    
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
    
    % middle node velocity along y axis
    v_m = zeros(num,1);
    midN = (N+1)/2;     % middle node number
    v_m(1) = u(2*midN); % update initial velocity
    
    for i=2:num
        % initial trial
        q=q0;
    
        % initialize error value
        err = 10*e;
    
        while err>e
            f=M/dt*((q-q0)/dt-u);  % discrete equation of motion
            J=M/dt^2;              % corresponding Jacobian
    
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
    
            % viscous force calculation
            % discrete viscous force
            f = f+C*(q-q0)/dt;
            J = J+C/dt; % corresponding Jacobian
    
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
        
        % update middle node velocity
        v_m(i) = u(2*midN);
    end

    % output: terminal velocity    
    n1 = v_m(num);
end