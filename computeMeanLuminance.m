
% Inputs:   1 or 2 picture matrixes read with imread. 
% Outputs:  luminances for 1 or 2 pictures, p value of the t-test testing the difference in luminance between two pictures if you have two 
function [luminance1, luminance2, p] = computeMeanLuminance(pic1, pic2, disps)
    gray1 = rgb2gray(pic1);
    luminance1 = mean(gray1(:));
    if ~isempty(pic2)
        gray2 = rgb2gray(pic2);
        luminance2 = mean(gray2(:));

        gray1_db = im2double(gray1);
        gray2_db = im2double(gray2);

        [~, p] = ttest2(gray1_db(:), gray2_db(:));
        % fprintf("luminance 1: %.2f - luminance 2: %.2f - p: %.4f\n", luminance1, luminance2, p)
        if disps; fprintf([inputname(1) ': %.2f - ' inputname(2) ': %.2f - p: %.4f\n'], luminance1, luminance2, p); end
    else
        luminance2 = [];
        p = [];
        if disps
            fprintf("luminance: %.2f\n", luminance1);
            fprintf([inputname(1) ': %.2f\n'], luminance1);
        end
    end
end