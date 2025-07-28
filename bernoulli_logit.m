function p =bernoulli_logit(x, a)
p = 1 ./ (1 + exp(-x .* (2 .* a - 1)));
end





