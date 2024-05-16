function cohens_d = cohend(x, y)

    mean_diff = mean(x) - mean(y);
    pooled_std = sqrt((std(x)^2 + std(y)^2) / 2);
    cohens_d = mean_diff / pooled_std;

end