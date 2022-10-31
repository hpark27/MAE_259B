% function to compute time parallel
% input: old reference frame, position vector q0, q
% output: reference frame ref1, ref2
function [ref1, ref2] = computeTimeParallel(ref1_old, q0, q)

% calculate number of edge
edge = (length(q)+1)/4-1;

% calculate tangent
tangent0 = computeTangent(q0); % old one : number of edge x 3
tangent = computeTangent(q);   % new one : number of edge x 3

% create reference frame director
ref1 = zeros(edge, 3); % first one
ref2 = zeros(edge, 3); % second one

% loop over edge
for i=1:edge
    t0 = tangent0(i,:); % old tangent on i-th edge
    t = tangent(i,:);   % new tangent on i-th edge

    ref1_local = parallel_transport(ref1_old(i,:), t0, t);

    % ref1 does not have a component parallel to tangent
    % enforce ref1 is perpendicular to t
    ref1_local = ref1_local - dot(ref1_local,t)*t;

    % enforce unit vector
    ref1_local = ref1_local / norm(ref1_local);

    % store the data
    ref1(i,:) = ref1_local;
    ref2(i,:) = cross(t, ref1_local); % cross product
end
end