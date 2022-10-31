% function to compute reference twist
% input: reference frame, tangent, reference twist
% output: rotated reference frame
function refTwist = computeRefTwist(ref1, tangent, refTwist)

% initialization
[edge,~] = size(ref1);

for i=2:edge % internal nodes
    u0 = ref1(i-1,:);    % ref1 vector of previous edge
    u1 = ref1(i,:);      % ref1 vector of current edge
    t0 = tangent(i-1,:); % tangent of previous edge
    t1 = tangent(i,:);   % tangent of current edge
    ut = parallel_transport(u0, t0, t1); % ut is different from u1

    ut = rotateAxisAngle(ut, t1, refTwist(i));
    refTwist(i) = refTwist(i) + signedAngle(ut, u1, t1);
end
end

% function to rotate axis angle
% input: vector v, axis z, angle theta
% output: rotated vector vNew
function vNew = rotateAxisAngle(v, z, theta)
if(theta==0)
    vNew = v;
else
    c = cos(theta);
    s = sin(theta);
    vNew = c*v + s*cross(z,v) + dot(z,v)*(1-c)*z;
end
end