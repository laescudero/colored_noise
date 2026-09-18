function result_matrix = octave_worker(in_data)
    % octave_worker  Helper function to run a single STFT realization in Octave parallel mode.
    %
    %   Usage:  result_matrix = octave_worker(in_data)
    %
    %   Input:
    %   in_data         :   a structure containing all necessary simulation parameters (e.g., L_sim, L, delta, sim_index).
    %
    %   Output:
    %   result_matrix   :   a 2-column matrix containing the detected zeros' coordinates [tau, omega].
    %
    % ---------------------------------------------------------
    addpath(genpath(fullfile(in_data.projectdir, 'includes')));

    % Rebuild the struct of the experiment
    experimentsParam = struct();
    experimentsParam.exp_id = in_data.exp_id;
    experimentsParam.use_mean = in_data.use_mean;
    experimentsParam.SNR = in_data.SNR;
    experimentsParam.saved_scale_factor = in_data.saved_scale_factor;

    % --- SEED & INDEX PASSING ---
    experimentsParam.sim_index = in_data.sim_index;
    if isfield(in_data, 'number_simulations')
        experimentsParam.number_simulations = in_data.number_simulations;
    end
    if isfield(in_data, 'use_fixed_seed')
        experimentsParam.use_fixed_seed = in_data.use_fixed_seed;
    end

    % Rebuild functions
    experimentsParam = build_experiment_functions(experimentsParam);

    [~, ~, ~, tau_zeros, omega_zeros] = run_single_realization(in_data.delta, in_data.L_sim, in_data.L, experimentsParam);

    result_matrix = [tau_zeros(:), omega_zeros(:)];
end
