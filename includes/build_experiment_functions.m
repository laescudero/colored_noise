function ep = build_experiment_functions(ep)
    % build_experiment_functions  Populates the experiment parameters struct with PSD and signal functions.
    %
    %   Usage:  ep = build_experiment_functions(ep)
    %
    %   Input:
    %   ep               :  a structure containing the experiment parameters.
    %
    %   Output:
    %   ep               :  the modified structure with added fields for PSD and signal.
    %
    %   Notation follows the manuscript "Zeros of the Spectrogram of Colored Noise":
    %   the parametric PSD model is S_{p,c}(omega) = K_{p,c} / (c + omega^2)^p (Eq. (5.26)),
    %   with K_{p,c} > 0 chosen so that the smoothed PSD satisfies (S_{p,c})_phi(0) = 1.
    %   The smoothing window is, as in the manuscript, phi(t) = e^{-2 pi t^2}.
    % ---------------------------------------------------------
    % Define the unnormalized PSD 1/(c+omega^2)^p and the deterministic signal f_1
    switch ep.exp_id
        case 1 % f_1 = 0, PSD S_{p,c} with p = 4, c = 0.5
            ep.PSD_unnormalized = @(t) (0.5 + t.^2).^(-4);
            ep.f1_func = @(t) 0;
        case 2 % f_1 = chirp, PSD S_{p,c} with p = 1, c = 1
            ep.PSD_unnormalized = @(t) (1 + t.^2).^(-1);
            ep.f1_func = @(t) exp(1i * pi * t.^2);
    end

    % Smoothing integral of the unnormalized PSD at omega = 0, k = 0.
    % This is (S_unnormalized)_phi (0); the normalization constant is K_{p,c} = 1/normalization_M0.
    M_k_aux = @(k, omega) arrayfun(@(w) integral(@(t) t.^k .* exp(-2 * pi * (w - t).^2) .* abs(ep.PSD_unnormalized(t)), w - 20, w + 20, 'RelTol', 1e-10, 'AbsTol', 1e-14), omega);
    normalization_M0 = M_k_aux(0, 0);

    % True normalized PSD S(omega), so that S_phi(0) = 1
    ep.PSD_S = @(t) (1 / normalization_M0) .* abs(ep.PSD_unnormalized(t));

    % Coloring filter transfer function sqrt(S(omega))
    ep.sqrt_PSD_S = @(t) sqrt(ep.PSD_S(t));

    % Moments of the smoothed PSD: M_k(omega) = \int t^k e^{-2 pi (omega - t)^2} S(t) dt.
    ep.M_k = @(k, omega) arrayfun(@(w) integral(@(t) t.^k .* exp(-2 * pi * (w - t).^2) .* ep.PSD_S(t), w - 20, w + 20, 'RelTol', 1e-12, 'AbsTol', 1e-15), omega);
end
