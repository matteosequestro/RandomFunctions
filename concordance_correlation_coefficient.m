function ccc = concordance_correlation_coefficient(x, y)
    mx = mean(x);
    my = mean(y);
    sx = std(x);
    sy = std(y);
    r = corr(x, y);

    ccc = (2 * r * sx * sy) / (sx^2 + sy^2 + (mx - my)^2);
end