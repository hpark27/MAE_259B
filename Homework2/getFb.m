% function to get bending energy
function [Fb, Jb] = getFb(q, m1, m2)
global kappaBar;
global EI;      % bending stiffness
global vorL;    % voroni length

node = (length(q)+1)/4;          % number of verticies
Fb = zeros(size(q));             % Force
Jb = zeros(length(q),length(q)); % Jacobian

for i=2:node-1
    node0 = [q(4*i-7), q(4*i-6), q(4*i-5)]; 
    node1 = [q(4*i-3), q(4*i-2), q(4*i-1)]; 
    node2 = [q(4*i+1), q(4*i+2), q(4*i+3)]; 

    m1e = m1(i-1,:); % first material director on c-1th edge
    m2e = m2(i-1,:); % second material director on c-1th edge
    m1f = m1(i,:);   % first material frame on cth edge
    m2f = m2(i,:);   % second material frame on cth edge

    [dF, dJ] = ...
        gradEb_hessEb(node0, node1, node2, ... 
        m1e, m2e, m1f, m2f, ...
        kappaBar(i,:), vorL(i), EI);

    index = 4*i-7:4*i+3; % index vector - size 11 

    % update force
    Fb(index) = Fb(index) - dF;

    % update Jacobian
    Jb(index, index) = Jb(index, index) - dJ;
end
end