
function [outs] = boxplotouts(x, criterion)
if ~exist('criterion', 'var'); criterion = 1.5; end
  q1 = quantile(x, 0.25);
  q3 = quantile(x, .75);
  iqr = q3 - q1;
  lb = q1 - (criterion * iqr);
  ub = q3 + (criterion * iqr);
  outs = x(x < lb | x > ub);
end

