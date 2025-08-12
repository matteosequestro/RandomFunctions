function p = bernoulli_logit(x)
    p = 1 ./ (1 + exp(-min(max(x,-50),50)));
end
