

% return the max number of identical consecutive elements (not sure the
% english is correct)
function [max_length] = samecons(x)

if width(x) == 1
    x = x';
end

diff_x = [1, diff(x)];
change_indices = find(diff_x ~= 0);
max_length = max(diff([0, change_indices, length(x) + 1]));


end