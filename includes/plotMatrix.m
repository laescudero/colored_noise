function plotMatrix(matrix, tau_coords, omega_coords, bartf)
    % plotMatrix  Log-scaled plot of a matrix with positive real coefficients.
    %   Usage:  plotMatrix(matrix, tau_coords, omega_coords, bartf);
    %           plotMatrix(matrix, tau_coords, omega_coords);
    %
    %   Input:
    %   matrix           :  a 2D matrix to be plotted (e.g. the spectrogram |V f|^2).
    %   tau_coords       :  an array with the time coordinates tau.
    %   omega_coords     :  an array with the frequency coordinates omega.
    %   bartf            :  a boolean indicating whether to plot a colorbar or not (default value: true).
    %   cRange           :  a variable that defines the range for the colorbar, (follows the convention [minRange, maxRange]).
    %
    %   Output:
    %   None.
    %
    % ---------------------------------------------------------
    global cRange

    if ~exist('bartf', 'var')
        bartf = true;
    end
    imagesc([tau_coords(1), tau_coords(length(tau_coords))], [omega_coords(length(omega_coords)), omega_coords(1)], rot90(20 * log10(matrix)), cRange);
    colormap(ltfat_inferno);
    axis xy;
    if bartf
        colorbar;
    end
    if ~bartf
        axis equal;
    end
    xlabel('Time ($\tau$)', 'Interpreter', 'latex');
    ylabel('Frequency ($\omega$)', 'Interpreter', 'latex');
    drawnow;
