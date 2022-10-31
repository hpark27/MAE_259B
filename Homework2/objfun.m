function [q, u, ref1, ref2] = objfun(q0, u, ref1, ref2, free, tol, ...
    refTwist)

% global variable
global Fg;   % gravitationanl force
global mMat; % mass matrix
global dt;   % time stamp

% initial guess
q = q0;

iteration = 1; % number of iteration

% error
err = 10*tol;

while err>tol
    % compute reference frame

    % temporary guess calculation
    [ref1_iterate, ref2_iterate] = computeTimeParallel(ref1,q0,q);
    
    % compute reference twist
    tangent = computeTangent(q); % compute tangent
    refTwist_iterate = computeRefTwist(ref1, tangent, refTwist);

    % compute material frame
    theta = q(4:4:end);
    [m1,m2] = computeMaterialDirectors(ref1_iterate, ref2_iterate, theta);

    % force and Jacobian calculation - Newton Raphson
    [Fb, Jb] = getFb(q, m1, m2);
    [Ft, Jt] = getFt(q, refTwist_iterate);
    [Fs, Js] = getFs(q);

    % force
    F = Fb + Ft + Fs + Fg;

    % jacobian of force
    jF = Jb + Jt + Js;
    
    % equations of motion
    f = mMat/dt*((q-q0)/dt-u)-F;

    % corresponding Jacobian
    J = mMat/dt^2-jF;

    % update free elements
    f_free = f(free);
    J_free = J(free, free);

    % Newton Raphson
    dq_free = J_free \ f_free;
    q(free) = q(free) - dq_free;

    % update error
    err = sum(abs(f_free));

    % update number of iteration
    iteration = iteration + 1;
end

% update velocity
u = (q - q0) / dt;

% update reference frame
ref1 = ref1_iterate;
ref2 = ref2_iterate;
end