% function to get twist energy
function [Ft, Jt] = getFt(q, refTwist)
global GJ;      % twisting stiffness
global vorL;    % voroni length

node = (length(q)+1)/4;          % number of verticies
Ft = zeros(size(q));             % Force
Jt = zeros(length(q),length(q)); % Jacobian

for i=2:node-1
    node0 = [q(4*i-7), q(4*i-6), q(4*i-5)]; 
    node1 = [q(4*i-3), q(4*i-2), q(4*i-1)]; 
    node2 = [q(4*i+1), q(4*i+2), q(4*i+3)]; 

    theta_e = q(4*i-4); % twist angle at i-1th edge
    theta_f = q(4*i);   % final angle at ith edge

    [dF, dJ] = ...
        gradEt_hessEt(node0, node1, node2, ... 
        theta_e, theta_f, refTwist(i),vorL(i), GJ);

    index = 4*i-7:4*i+3; % index vector - size 11 

    % update force
    Ft(index) = Ft(index) - dF;

    % update Jacobian
    Jt(index, index) = Jt(index, index) - dJ;
end
end