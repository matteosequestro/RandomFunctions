function out = movmean_upto(x, winlen)

lenx = length(x);
out = ones(1, lenx );
for ii = 1 : lenx 
    
    if ii > winlen+1
        out(ii) = mean(x(ii-winlen : ii));
    else
        out(ii) = mean(x(1 : ii));
    end
    
end

end