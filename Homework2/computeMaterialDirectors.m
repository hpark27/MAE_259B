% function to compute material director
% input: reference frame1, reference frame2, angle
% output: material director m1 and m2
function [m1,m2] = computeMaterialDirectors(ref1, ref2, theta)
% number of edge
edge = length(theta);

% material frame directors for all edge
m1 = zeros(edge,3); % first one
m2 = zeros(edge,3); % second one

% update frame director elements. Loop over edges
for i=1:edge
    c = cos(theta(i)); % cosine
    s = sin(theta(i)); % sine

    % update elements
    m1(i,:) = c*ref1(i,:)+s*ref2(i,:);
    m2(i,:) = -s*ref1(i,:)+c*ref2(i,:);
end
end