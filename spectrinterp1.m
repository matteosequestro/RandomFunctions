% based on: Mewett, D.T., Reynolds, K.J., Nazeran, H., 2004. Reducing power line interference in
% digitised electromyogram recordings by spectrum interpolation. Medical and
% Biological Engineering and Computing 42 (4), 524–531. https://doi.org/10.1007/
% BF02350994.

function [cleansignal, Y] = spectrinterp1(signal, frequencies, fs)
    
    % Do fft
    Y = fft(signal);

    % You will add random noise to the interpolated values (for the moment)
    noise = 0.2;    

    % Loop through all the frequencies to interpolate
    for i = 1 : height(frequencies)

        % Find the position of the frequencies on the FFT
        toInterValues = round(frequencies(i, :) .* length(Y) / fs, 0);
        
        % Interpolate values for left side
        interp1 = linspace(Y(toInterValues(1) -1), ...
                            Y(toInterValues(2) + 1), ...
                            length( Y(toInterValues(1) : toInterValues(2)) ));
        
        % Interpolate values for right side
        interp2 = linspace(Y(end - (toInterValues(1) -1)), ...
                            Y(end - (toInterValues(2) + 1)), ...
                            length( Y(end - (toInterValues(1) : toInterValues(2))) ));
        
        
        % Replace interpolated values on left and right sides of fft
        Y(toInterValues(1) : toInterValues(2)) = interp1;
        Y(end - (toInterValues(1) : toInterValues(2))) = interp2;
    end

    % Reconstruct the signal with inverse fft
    cleansignal = ifft(Y,'symmetric');

end