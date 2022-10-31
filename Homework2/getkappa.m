% function to get kappa
% input: position vector q, material director m1, m2
% output: kappa k
function k = getkappa(q, m1, m2)
node = (length(q)+1)/4; % same with computeTangent function
edge = node-1;          % same with computeTangent function

k = zeros(node,2); % initialize kappa

% start from the second node
% first and second node do not have curvature
for i=2:edge
    node0 = q(4*i-7:4*i-5); 
    node1 = q(4*i-3:4*i-1); 
    node2 = q(4*i+1:4*i+3);

    m1e = m1(i-1, :); % m1 vector on (i-1)th edge
    m2e = m2(i-1, :); % m2 vector on (i-1)th edge
    m1f = m1(i, :); % m1 vector on ith edge
    m2f = m2(i, :); % m2 vector on ith edge

    % get the curvature for each node
    kLocal = computekappa(node0,node1, node2, m1e, m2e, m1f, m2f );

    % store the derived value into matrix
    k(i,1) = kLocal(1);
    k(i,2) = kLocal(2);
end
end