function [sw] = slidingWindowBinary(x,y,windowWidth,timeMax)

   for i = 1:timeMax

        lower = i - (windowWidth/2);
        upper = i + (windowWidth/2);

        %processing for unchanged trials
        points = find(x*1000 <=upper & x*1000 >=lower);
        sample = y(points,:);

        if isempty(sample) == 0 
        %set sliding window results:
            sw(i) = sum(sample)/length(sample); 
        else
            sw(i) = NaN; 
        end
   end
end