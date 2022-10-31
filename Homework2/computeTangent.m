% Function to compute tangent
% Input parameter : a vector with size 4*node-1
% Output : Tangent in the form of a matrix with size (edge, 3)
% where edge = node-1
function tan = computeTangent(q)
node = (length(q)+1)/4;
edge = node-1;

% initialize matrix
tan = zeros(node,3);

% update matrix elements
for i=1:edge % loop over edge
    xc = q(4*i-3:4*i-1);
    xcp1 = q(4*i+1:4*i+3);
    dx=xcp1-xc;
    tan(i,:)=dx/norm(dx);
end
end