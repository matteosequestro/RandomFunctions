% maximum number of time a number (crit) is repeated in consecutive 
% positions in a vector (bad english I know I'm tired)
function maxcon = maxcons(x, crit)

wherecon = find(x == crit);
wherecon2 = diff(wherecon);
ncon = 0; maxcon = 0;
for ii = 1  :length(wherecon2)
    if wherecon2(ii) == 1
        ncon = ncon +1;
    else
        ncon = 0;
    end
    if ncon > maxcon; maxcon = ncon; end
end



end