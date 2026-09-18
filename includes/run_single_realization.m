function [STFT_analysis, tau_analysis, omega_analysis, tau_zeros, omega_zeros] = run_single_realization(delta, L_sim, L, experimentsParam)
    % run_single_realization  Computes a single STFT realization of colored noise (and signal).
    %
    %   Usage:  [STFT_analysis, tau_analysis, omega_analysis, tau_zeros, omega_zeros] = run_single_realization(delta, L_sim, L, experimentsParam)
    %
    %   Input:
    %   delta            :  grid step size delta = 1/(2*L_sim) for the STFT.
    %   L_sim            :  simulation half-width (STFT computed over [-L_sim,L_sim)).
    %   L                :  analysis half-width (zeros/STFT restricted to [-L,L]^2).
    %   experimentsParam :  struct containing PSD and deterministic signal functions.
    %
    %   Output:
    %   STFT_analysis    :  the computed STFT V(f) restricted to the analysis grid.
    %   tau_analysis     :  time coordinates corresponding to the STFT matrix.
    %   omega_analysis   :  frequency coordinates corresponding to the STFT matrix.
    %   tau_zeros        :  time coordinates of the detected zeros.
    %   omega_zeros      :  frequency coordinates of the detected zeros.
    %
    % ---------------------------------------------------------

    tau_grid = -L_sim:delta:(L_sim - delta);
    fs = 1 / delta;

    % --- REPRODUCIBILITY (SEED CONTROL) ---
    if isfield(experimentsParam, 'sim_index')
        if isfield(experimentsParam, 'use_fixed_seed') && experimentsParam.use_fixed_seed
            % FIXED SEED: Seeds like in the paper
            n_val = 0;
            if isfield(experimentsParam, 'number_simulations')
                n_val = experimentsParam.number_simulations;
            end

            seed_val = n_val * 100000 + experimentsParam.sim_index * 1000;

            if exist('OCTAVE_VERSION', 'builtin') ~= 0
                rand('seed', seed_val);
                randn('seed', seed_val);
            else
                rng(seed_val, 'twister');
            end
        elseif exist('OCTAVE_VERSION', 'builtin') ~= 0
            % OCTAVE RANDOM SEED: Ensure completely different seeds across processes
            seed_val = experimentsParam.sim_index * 1000 + mod(round(cputime * 10000), 10000);
            rand('seed', seed_val);
            randn('seed', seed_val);
        end
    end

    % Generate complex white noise zeta = (xi_1 + i xi_2)/sqrt(2)
    white_noise = 1 / sqrt(2) .* (randn(size(tau_grid)) + 1i .* randn(size(tau_grid)));
    freq_grid = linspace(-fs / 2, fs / 2 - delta, length(white_noise));

    % Color the noise with the filter sqrt(S(omega)) to obtain the discrete colored noise W
    sqrt_PSD_vals = ifftshift(experimentsParam.sqrt_PSD_S(freq_grid));
    colored_noise = ifft(fft(white_noise) .* sqrt_PSD_vals);

    % Add the deterministic part f_1 using the scale factor : f = f_1 + W
    if isfield(experimentsParam, 'use_mean') && experimentsParam.use_mean
        f1_vals = experimentsParam.f1_func(tau_grid);
        signal = experimentsParam.saved_scale_factor .* f1_vals + colored_noise;
    else
        signal = colored_noise;
    end

    [discreteSTFT, tau_grid, omega_grid] = paper_stft(signal, L_sim, delta);

    % Boundary cropping with a buffer of 1 bin.
    margin_phys = delta;
    valid_tau_pad = (tau_grid >= -L - margin_phys) & (tau_grid <= L + margin_phys);
    valid_omega_pad = (omega_grid >= -L - margin_phys) & (omega_grid <= L + margin_phys);

    STFT_padded = discreteSTFT(valid_omega_pad, valid_tau_pad);
    tau_pad = tau_grid(valid_tau_pad);
    omega_pad = omega_grid(valid_omega_pad);

    [locminomegaind, locmintauind] = MGN(abs(STFT_padded).^2);
    tau_zeros = tau_pad(locmintauind);
    omega_zeros = omega_pad(locminomegaind);

    val_tau_analysis = (tau_grid >= -L) & (tau_grid <= L);
    val_omega_analysis = (omega_grid >= -L) & (omega_grid <= L);
    STFT_analysis = discreteSTFT(val_omega_analysis, val_tau_analysis);
    tau_analysis = tau_grid(val_tau_analysis);
    omega_analysis = omega_grid(val_omega_analysis);
end
