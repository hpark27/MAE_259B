% function to get stretch energy
function [Fs, Js] = getFs(q)
global EA;      % stretch stifness
global refL;    % reference length

node = (length(q)+1)/4;          % number of verticies
edge = node - 1;                 % edge
Fs = zeros(size(q));             % Force
Js = zeros(length(q),length(q)); % Jacobian

for i=1:edge
    node0 = [q(4*i-3), q(4*i-2), q(4*i-1)]; 
    node1 = [q(4*i+1), q(4*i+2), q(4*i+3)]; 

    [dF, dJ] = gradEs_hessEs(node0, node1, refL(i),EA);

     % index vector - size 6
    index = [4*i-3, 4*i-2, 4*i-1, 4*i+1, 4*i+2, 4*i+3];

    % update force
    Fs(index) = Fs(index) - dF;

    % update Jacobian
    Js(index, index) = Js(index, index) - dJ;
end
end