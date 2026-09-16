function demo_spectrogram()
    % demo_spectrogram  Computes and plots the spectrogram and its zeros for a signal with colored noise.
    %   Also computes, for comparison: (i) the spectrogram of the same
    %   deterministic signal with the original (uncolored) white noise, and
    %   (ii) the PSDs of both the colored and the white noise.
    %
    %   Usage:  demo_spectrogram()
    %
    %   Input:
    %   None. Parameters are defined inside the function. See Section "1. User Parameters" below.
    %
    %   Output:
    %   None. Generates and saves a plot in the 'images' folder.
    %
    % ---------------------------------------------------------
    clear all;
    close all;

    % Ensure auxiliary functions are in the path
    addpath(genpath(fullfile(pwd, 'includes')));

    % =========================================================================
    % 1. USER PARAMETERS (The user is invited to modify these parameters)
    % =========================================================================
    SNR_dB = 2; % Signal-to-Noise Ratio (SNR) in dB

    % Deterministic signal f_1(t)
    f1_func = @(t) exp(-pi * t.^2);

    % Power Spectral Density of the noise S(\omega)
    PSD_S = @(w) (0.01 + 2 * (w - 2.5).^4).^(-1);

    % Output configuration
    output_folder = fullfile(pwd, 'images');
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % =========================================================================
    % 2. GRID CONFIGURATION
    % =========================================================================
    L_sim = 40;             % STFT calculation support.
    L = 20;                 % Math bounding box for zeros and spectrogram.
    L_plot = 6;             % Visual bounding box for the final plot to zoom-in.

    delta = 1 / (2 * L_sim);
    fs = 1 / delta;

    tau_grid = -L_sim:delta:(L_sim - delta);

    disp('Generating signals and calculating spectrogram...');

    % =========================================================================
    % 3. NOISE GENERATION AND SNR SCALING
    % =========================================================================
    % Evaluate the square root of the PSD for the filter
    freq_grid = linspace(-fs / 2, fs / 2 - delta, length(tau_grid));
    sqrt_PSD_vals = ifftshift(sqrt(PSD_S(freq_grid)));

    % Expected squared L^2-Norm of the noise
    expected_power_noise = mean(abs(sqrt_PSD_vals).^2);

    % Evaluate squared L^2-Norm of the deterministic component
    f1_vals = f1_func(tau_grid);
    power_sig = mean(abs(f1_vals).^2);

    % Scaling to achieve SNR
    if power_sig > 0
        scale_factor = sqrt(10^(SNR_dB / 10) * expected_power_noise / power_sig);
    else
        scale_factor = 1;
    end

    % Generate a realization of white noise and color it
    white_noise = 1 / sqrt(2) .* (randn(size(tau_grid)) + 1i .* randn(size(tau_grid)));
    colored_noise = ifft(fft(white_noise) .* sqrt_PSD_vals);

    % Final signal = Scaled deterministic signal + Colored noise
    final_signal = (scale_factor .* f1_vals) + colored_noise;

    % Same deterministic signal, but with the original (uncolored) white noise
    signal_white_noise = (scale_factor .* f1_vals) + white_noise;

    % =========================================================================
    % 4. STFT AND ZERO ESTIMATION
    % =========================================================================
    [STFT_mat, tau_grid, omega_grid] = paper_stft(final_signal, L_sim, delta);
    Spectrogram_mag = abs(STFT_mat).^2;

    % Crop of the STFT domain to [-L, L]^2 with a margin of a bin.
    margin = delta;
    val_tau_pad = (tau_grid >= -L - margin) & (tau_grid <= L + margin);
    val_omega_pad = (omega_grid >= -L - margin) & (omega_grid <= L + margin);

    STFT_padded = Spectrogram_mag(val_omega_pad, val_tau_pad);
    tau_pad = tau_grid(val_tau_pad);
    omega_pad = omega_grid(val_omega_pad);

    % Find local minima (zeros)
    [locminomegaind, locmintauind] = MGN(STFT_padded);
    zeros_tau = tau_pad(locmintauind);
    zeros_omega = omega_pad(locminomegaind);

    % Crop grids for the final plot data mapping
    val_tau_plot = (tau_grid >= -L) & (tau_grid <= L);
    val_omega_plot = (omega_grid >= -L) & (omega_grid <= L);
    Spectrogram_plot = Spectrogram_mag(val_omega_plot, val_tau_plot);
    tau_plot = tau_grid(val_tau_plot);
    omega_plot = omega_grid(val_omega_plot);

    % =========================================================================
    % 4bis. SPECTROGRAM OF THE SIGNAL WITH THE ORIGINAL WHITE NOISE
    % =========================================================================
    [STFT_mat_wn, ~, ~] = paper_stft(signal_white_noise, L_sim, delta);
    Spectrogram_mag_wn = abs(STFT_mat_wn).^2;

    STFT_padded_wn = Spectrogram_mag_wn(val_omega_pad, val_tau_pad);

    % Find local minima (zeros)
    [locminomegaind_wn, locmintauind_wn] = MGN(STFT_padded_wn);
    zeros_tau_wn = tau_pad(locmintauind_wn);
    zeros_omega_wn = omega_pad(locminomegaind_wn);

    % Same time/frequency crop as the main spectrogram, so both panels share axes.
    Spectrogram_wn_plot = Spectrogram_mag_wn(val_omega_plot, val_tau_plot);

    % =========================================================================
    % 4ter. PSD OF THE NOISE
    % =========================================================================
    % Normalized by its own maximum, purely for display, so it shares the
    % same [0, 1] scale as the white noise's PSD below (the un-normalized
    % PSD_S is what actually defines the coloring filter sqrt_PSD_vals above).
    psd_vals = PSD_S(omega_plot);
    psd_vals_norm = psd_vals / max(psd_vals);

    % PSD of the original (uncolored) white noise: flat, S(omega) = 1.
    psd_vals_wn = ones(size(omega_plot));

    % =========================================================================
    % 5. VISUALIZATION AND EXPORT
    % =========================================================================
    % 1. cRange
    global cRange
    cRange = [-100, -10];

    AXES_POS = [0.15, 0.15, 0.68, 0.75];

    % --- Figure 1: Spectrogram of the signal + colored noise, with its zeros ---
    fig1 = figure('Name', 'Demo: Spectrogram and Zeros of Colored Noise', 'Position', [100, 100, 800, 600]);
    plotMatrix(Spectrogram_plot.', tau_plot, omega_plot);
    set(gca, 'YDir', 'normal');
    colorbar;
    hold on;
    scatter(zeros_tau, zeros_omega, 15, 'w', 'filled', 'MarkerEdgeColor', 'k');
    hold off;
    xlim([-L_plot, L_plot]);
    ylim([-L_plot, L_plot]);
    xlabel('Time (\tau)', 'Interpreter', 'tex', 'FontSize', 16);
    ylabel('Frequency (\omega)', 'Interpreter', 'tex', 'FontSize', 16);
    set(gca, 'Position', AXES_POS);

    output_filename = fullfile(output_folder, sprintf('demo_spectrogram.png'));
    export_figure(fig1, output_filename);
    disp(['Figure successfully saved to: ', output_filename]);

    % --- Figure 2: Spectrogram of the signal with the original (uncolored) white noise ---
    fig2 = figure('Name', 'Demo: Spectrogram with the Original White Noise', 'Position', [100, 100, 800, 600]);
    plotMatrix(Spectrogram_wn_plot.', tau_plot, omega_plot);
    set(gca, 'YDir', 'normal');
    colorbar;
    hold on;
    scatter(zeros_tau_wn, zeros_omega_wn, 15, 'w', 'filled', 'MarkerEdgeColor', 'k');
    hold off;
    xlim([-L_plot, L_plot]);
    ylim([-L_plot, L_plot]);
    xlabel('Time (\tau)', 'Interpreter', 'tex', 'FontSize', 16);
    ylabel('Frequency (\omega)', 'Interpreter', 'tex', 'FontSize', 16);
    set(gca, 'Position', AXES_POS);

    output_filename_wn = fullfile(output_folder, sprintf('demo_spectrogram_whitenoise.png'));
    export_figure(fig2, output_filename_wn);
    disp(['Figure successfully saved to: ', output_filename_wn]);

    % --- Figure 3: PSD of the noise, S(omega) vs omega (normalized by its max) ---
    fig3 = figure('Name', 'Demo: PSD of the Noise', 'Position', [100, 100, 800, 600]);
    plot(omega_plot, psd_vals_norm, 'LineWidth', 1.5);
    grid on;
    xlim([-L_plot, L_plot]);
    ylim([0, 1.2]);
    xlabel('Frequency (\omega)', 'Interpreter', 'tex', 'FontSize', 16);
    ylabel('S(\omega) / max S', 'Interpreter', 'tex', 'FontSize', 16);
    set(gca, 'Position', AXES_POS);

    output_filename_psd = fullfile(output_folder, sprintf('demo_spectrogram_psd.png'));
    export_figure(fig3, output_filename_psd);
    disp(['Figure successfully saved to: ', output_filename_psd]);

    % --- Figure 4: PSD of the white noise, S(omega) = 1 vs omega ---
    fig4 = figure('Name', 'Demo: PSD of the White Noise', 'Position', [100, 100, 800, 600]);
    plot(omega_plot, psd_vals_wn, 'LineWidth', 1.5);
    grid on;
    xlim([-L_plot, L_plot]);
    ylim([0, 1.2]);
    xlabel('Frequency (\omega)', 'Interpreter', 'tex', 'FontSize', 16);
    ylabel('S(\omega) / max S', 'Interpreter', 'tex', 'FontSize', 16);
    set(gca, 'Position', AXES_POS);

    output_filename_psd_wn = fullfile(output_folder, sprintf('demo_spectrogram_psd_whitenoise.png'));
    export_figure(fig4, output_filename_psd_wn);
    disp(['Figure successfully saved to: ', output_filename_psd_wn]);
end

function export_figure(fig, filename)
    drawnow;
    frame = getframe(fig);
    img = frame.cdata;

    border_color = double(reshape(img(1, 1, :), 1, 3));
    if ~isequal(border_color, [255, 255, 255])
        mask = all(double(img) == reshape(border_color, 1, 1, 3), 3);
        img(repmat(mask, [1, 1, 3])) = 255;
    end

    imwrite(img, filename);
end
