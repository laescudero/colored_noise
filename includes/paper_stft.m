function [stft_result, tau, omega] = paper_stft(signal, L_sim, delta)
    % PAPER_STFT  Discrete STFT of Eq. (5.12), "Zeros of the Spectrogram of Colored Noise":
    %   V_phi f(tau_k,omega_j) = delta * sum_s 2^(1/4) f_s
    %                            * exp(-pi(tau_s-tau_k)^2) exp(-2pi i tau_s omega_j),
    %   with tau_s, tau_k, omega_j = -L_sim + (s,k,j)*delta.
    %
    %   The omega_j sum is evaluated with Bluestein's algorithm.
    %   See https://ccrma.stanford.edu/~jos/mdft/Bluestein_s_FFT_Algorithm.html
    %
    %   Input:
    %   signal : sampled signal f(s).
    %   L_sim  : half-width of the simulation time interval [-L_sim,L_sim)
    %   delta  : sampling step (delta = 1/(2*L_sim))
    %
    %   Output:
    %   stft_result(j,k) : discrete STFT V(f) as in Eq. (5.12).
    %   tau              : mapping from integer to time coordinates
    %   omega            : mapping from integer to frequency coordinates
    %
    % =========================================================================

    signal = signal(:).';
    K = length(signal);
    tau = -L_sim:delta:(L_sim - delta);
    omega = tau;

    stft_result = zeros(K, K);

    % =========================================================================
    % PRECOMPUTATIONS
    % =========================================================================
    s = 0:(K - 1);
    j = 0:(K - 1);

    alpha = delta^2;
    W_pre = exp(-1i * pi * alpha * s.^2);
    W_post = exp(-1i * pi * alpha * j.^2);
    W_conv = exp(1i * pi * alpha * ((-K + 1):(K - 1)).^2);

    % Pad to the next power of two for a faster FFT.
    N_fft = 2^nextpow2(2 * K - 1);
    fft_W_conv = fft(W_conv, N_fft);

    % Shift phases due to integration starting at -L_sim instead of 0
    phase_s = exp(2i * pi * s * delta * L_sim);
    phase_j = exp(2i * pi * j * delta * L_sim);
    phase_const = exp(-2i * pi * L_sim^2);

    % =========================================================================
    % MAIN LOOP
    % =========================================================================
    for k = 1:K
        center = tau(k);

        % 1. Apply normalized and time-shifted Gaussian window phi(t) = 2^(1/4) e^{-pi t^2}
        win = 2^(1 / 4) * exp(-pi * (tau - center).^2);
        v = signal .* win .* delta;

        % 2. Prepare signal for convolution
        X = v .* phase_s;
        A = X .* W_pre;

        % 3. FFT
        fft_A = fft(A, N_fft);
        conv_res = ifft(fft_A .* fft_W_conv);

        % 4. Extract valid convolution overlap
        S = conv_res(K:2 * K - 1);

        % 5. Correction of the phase obtain the discrete STFT as in Eq. (5.12)
        S = S .* W_post;
        stft_result(:, k) = (S .* phase_j .* phase_const).';
    end
end
