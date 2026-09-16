function main()
    % main  Entry point for running STFT zero simulations and plotting results.
    %      It computes the simulations described in Section 5 and 6 of the
    %      manuscript "Zeros of the Spectrogram of Colored Noise".
    %
    %   Usage:  main()
    %
    %   Input:
    %   None. Check the definition of build_experiment_functions if you want to
    %         generate other experiments.
    %         You can also set experimentsParam.use_fixed_seed = false; if you
    %         do not want to run the same seeds as we did for our manuscript.
    %
    %   Output:
    %   Depending on the chosen menu option, it generates:
    %     - Option 1: Runs simulations and saves the raw data to a .mat file
    %                 in the project root directory.
    %     - Option 2: Loads the .mat file, computes the PSD and density
    %                 estimators, generates the plots, and exports
    %                 the numerical results to .dat files inside the 'images' folder.
    %     - Option 3: Performs both Option 1 and Option 2 sequentially.
    %
    % ---------------------------------------------------------
    % --- Initial Parameters ---
    clear all;
    projectdir = pwd;
    addpath(genpath(fullfile(projectdir, 'includes')));

    folderPath = fullfile(projectdir, 'images');
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    disp(['Output folder: ', folderPath]);

    L = 25;
    L_sim = 50;
    delta = 1 / (2 * L_sim);

    fprintf('Resolution of the grid: %d\n', delta);
    fprintf('L: %d\n', L);

    array_number_simulations = [5, 10, 100, 1000];

    % --- Menu ---
    disp('Please choose the experiment to run/plot:');
    disp('1: Exp 1 - ZERO MEAN -  PSD: S_{p,c} with p=4, c=0.5');
    disp('2: Exp 2 - CHIRP     -  PSD: S_{p,c} with p=1, c=1');
    disp('3: RUN ALL EXPERIMENTS');
    experiment_choice = input('Enter your choice (1-3): ');

    if experiment_choice == 3
        experiments_to_run = 1:2;
    elseif ismember(experiment_choice, 1:2)
        experiments_to_run = experiment_choice;
    else
        disp('Invalid choice.');
        return
    end

    disp(' ');
    disp('Please choose an action:');
    disp('1: Run simulations and save raw zeros data');
    disp('2: Load data and generate plots and export .dat files');
    disp('3: Run simulations, save data, generate plots and export .dat files');
    choice = input('Enter your choice (1, 2, or 3): ');

    % =========================================================================
    % EXPERIMENT LOOP
    % =========================================================================
    for curr_exp = experiments_to_run

        % Seed Reset for reproducibility
        try
            rng("default");
        catch
            rand('seed', 0);
            randn('seed', 0);
        end

        experimentsParam = struct();
        experimentsParam.exp_id = curr_exp;
        experimentsParam.use_fixed_seed = true;

        if curr_exp == 1
            experimentsParam.use_mean = false;
            file_prefix = 'exp1';
        elseif curr_exp == 2
            experimentsParam.use_mean = true;
            experimentsParam.SNR = 20;
            file_prefix = 'exp2';
        end

        output_filename = fullfile(projectdir, sprintf('simulation_results_%s.mat', file_prefix));

        if choice == 1 || choice == 3
            run_simulations_and_save(projectdir, folderPath, L_sim, L, delta, ...
                                     experimentsParam, array_number_simulations, output_filename, file_prefix);
        end

        if choice == 2 || choice == 3
            if exist(output_filename, 'file')
                load_data_and_plot(folderPath, output_filename, file_prefix);
            else
                disp(['Error: ', output_filename, ' not found.']);
            end
        end
    end
end

% =========================================================================
% MODULE 1: SIMULATION AND SAVING
% =========================================================================
function run_simulations_and_save(projectdir, folderPath, L_sim, L, delta, ...
                                  experimentsParam, array_number_simulations, output_filename, file_prefix)
    % run_simulations_and_save  Executes the simulation loop for STFT zero counting and saves raw results.
    %
    %   Usage:  run_simulations_and_save(projectdir, folderPath, L_sim, L, delta, ...
    %                                    experimentsParam, array_number_simulations, output_filename, file_prefix)
    %
    %   Input:
    %   projectdir               :  base directory path of the project.
    %   folderPath               :  output folder path.
    %   L_sim                    :  simulation half-width (STFT computed over [-L_sim,L_sim)).
    %   L                        :  analysis half-width (zeros counted on [-L,L]^2).
    %   delta                    :  grid resolution for time and frequency (delta = 1/(2*L_sim)).
    %   experimentsParam         :  struct containing experimental parameters.
    %   array_number_simulations :  array indicating numbers of realizations n to run.
    %   output_filename          :  path to save the .mat file.
    %   file_prefix              :  prefix for file naming.
    %
    %   Output:
    %   None. Saves the variables to a MAT file.
    %
    % ---------------------------------------------------------

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if ~isOctave
        pool = gcp('nocreate');
        if isempty(pool)
            parpool();
        end
    end

    experimentsParam = build_experiment_functions(experimentsParam);

    tau_grid = -L_sim:delta:(L_sim - delta);
    fs = 1 / delta;

    % Precalculate scale factor to achieve the SNR
    if isfield(experimentsParam, 'use_mean') && experimentsParam.use_mean
        freq_grid_base = linspace(-fs / 2, fs / 2 - delta, length(tau_grid));
        sqrt_PSD_base = ifftshift(experimentsParam.sqrt_PSD_S(freq_grid_base));
        expected_power_noise = mean(abs(sqrt_PSD_base).^2);

        f1_vals = experimentsParam.f1_func(tau_grid);
        power_sig = mean(abs(f1_vals).^2);

        if power_sig > 0
            experimentsParam.saved_scale_factor = sqrt(10^(experimentsParam.SNR / 10) * expected_power_noise / power_sig);
        else
            experimentsParam.saved_scale_factor = 1;
        end
    end

    resultsArray = struct();

    progress_counter = 0;
    total_sims = 0;

    for sims_idx = 1:length(array_number_simulations)
        n_sims = array_number_simulations(sims_idx);
        fprintf('\n[%s] Experiment running with n = %d.\n', upper(strrep(file_prefix, '_', ' ')), n_sims);

        all_omega_zeros_cell = cell(n_sims, 1);
        all_tau_zeros_cell = cell(n_sims, 1);

        % ITERATION 1 (Computes the base realization)
        experimentsParam.sim_index = 1;
        experimentsParam.number_simulations = array_number_simulations(sims_idx);
        [stft_mat, tau_analysis, omega_analysis, tau_zeros, omega_zeros] = run_single_realization(delta, L_sim, L, experimentsParam);

        resultsArray(sims_idx).L = L;
        resultsArray(sims_idx).number_of_simulations = n_sims;

        if sims_idx == 1
            resultsArray(sims_idx).first_stft = stft_mat;
            resultsArray(sims_idx).first_tau_zeros = tau_zeros;
            resultsArray(sims_idx).first_omega_zeros = omega_zeros;
        else
            resultsArray(sims_idx).first_stft = [];
            resultsArray(sims_idx).first_tau_zeros = [];
            resultsArray(sims_idx).first_omega_zeros = [];
        end

        resultsArray(sims_idx).tau_grid = tau_analysis;
        resultsArray(sims_idx).omega_grid = omega_analysis;

        all_omega_zeros_cell{1} = omega_zeros(:);
        all_tau_zeros_cell{1} = tau_zeros(:);

        % ITERATIONS 2 to n
        if n_sims > 1
            if ~isOctave
                progressQueue = parallel.pool.DataQueue;
                afterEach(progressQueue, @(~) update_progress_display());
                progress_counter = 1;
                total_sims = n_sims;

                parfor index = 2:n_sims
                    ep_local = experimentsParam;
                    ep_local.sim_index = index;
                    [~, ~, ~, ltau, lomega] = run_single_realization(delta, L_sim, L, ep_local);
                    all_omega_zeros_cell{index} = lomega(:);
                    all_tau_zeros_cell{index} = ltau(:);
                    send(progressQueue, true);
                end
            else
                % --- OCTAVE PARALLEL EXECUTION ---
                try
                    pkg load parallel;
                    n_cores = nproc();
                    fprintf('Running Octave in PARALLEL mode using %d cores...\n', n_cores);

                    worker_inputs = cell(1, n_sims - 1);
                    for k = 1:(n_sims - 1)
                        sim_real_idx = k + 1;
                        in_data = struct();
                        in_data.sim_index = sim_real_idx;
                        in_data.projectdir = projectdir;
                        in_data.L_sim = L_sim;
                        in_data.L = L;
                        in_data.delta = delta;
                        in_data.exp_id = experimentsParam.exp_id;
                        in_data.number_simulations = n_sims;

                        if isfield(experimentsParam, 'use_fixed_seed')
                            in_data.use_fixed_seed = experimentsParam.use_fixed_seed;
                        end

                        if isfield(experimentsParam, 'use_mean') && experimentsParam.use_mean
                            in_data.use_mean = true;
                            in_data.SNR = experimentsParam.SNR;
                            if isfield(experimentsParam, 'saved_scale_factor')
                                in_data.saved_scale_factor = experimentsParam.saved_scale_factor;
                            else
                                in_data.saved_scale_factor = 1;
                            end
                        else
                            in_data.use_mean = false;
                            in_data.SNR = 0;
                            in_data.saved_scale_factor = 1;
                        end
                        worker_inputs{k} = in_data;
                    end

                    parallel_results = parcellfun(n_cores, "octave_worker", worker_inputs, "UniformOutput", false, "VerboseLevel", 1);

                    for i = 1:length(parallel_results)
                        real_idx = i + 1;
                        matrix_data = parallel_results{i};
                        all_tau_zeros_cell{real_idx}   = matrix_data(:, 1);
                        all_omega_zeros_cell{real_idx} = matrix_data(:, 2);
                    end

                catch err
                    disp(['Parallel error: ', err.message]);
                    disp('Running Octave in SERIAL mode...');
                    for index = 2:n_sims
                        experimentsParam.sim_index = index;
                        [~, ~, ~, ltau, lomega] = run_single_realization(delta, L_sim, L, experimentsParam);
                        all_omega_zeros_cell{index} = lomega(:);
                        all_tau_zeros_cell{index} = ltau(:);
                        fprintf('Progress: %d / %d\n', index, n_sims);
                    end
                end
            end
        end
        resultsArray(sims_idx).omega_zeros_all = vertcat(all_omega_zeros_cell{:});
        resultsArray(sims_idx).tau_zeros_all = vertcat(all_tau_zeros_cell{:});
        resultsArray(sims_idx).omega_zeros_by_realization = all_omega_zeros_cell;
    end

    if isfield(experimentsParam, 'PSD_unnormalized')
        experimentsParam = rmfield(experimentsParam, 'PSD_unnormalized');
    end
    if isfield(experimentsParam, 'f1_func')
        experimentsParam = rmfield(experimentsParam, 'f1_func');
    end
    if isfield(experimentsParam, 'sqrt_PSD_S')
        experimentsParam = rmfield(experimentsParam, 'sqrt_PSD_S');
    end
    if isfield(experimentsParam, 'PSD_S')
        experimentsParam = rmfield(experimentsParam, 'PSD_S');
    end
    if isfield(experimentsParam, 'M_k')
        experimentsParam = rmfield(experimentsParam, 'M_k');
    end

    if isOctave
        save(output_filename, 'resultsArray', 'experimentsParam', 'L_sim', 'L', 'delta', 'array_number_simulations', '-v7');
    else
        save(output_filename, 'resultsArray', 'experimentsParam', 'L_sim', 'L', 'delta', 'array_number_simulations', '-v7.3');
    end

    disp(['Data saved to: ', output_filename]);

    function update_progress_display()
        progress_counter = progress_counter + 1;
        if mod(progress_counter, 10) == 0
            fprintf('Progress: %d / %d complete.\n', progress_counter, total_sims);
        end
    end

end

% =========================================================================
% MODULE 2: LOADING, COMPUTING ESTIMATORS, PLOTTING, AND EXPORT (.dat)
% =========================================================================
function load_data_and_plot(folderPath, data_filename, file_prefix)
    % load_data_and_plot  Computes estimators, exports to .dat files, and plots results.
    %
    %   Usage:  load_data_and_plot(folderPath, data_filename, file_prefix)
    %
    %   Input:
    %   folderPath      :   output folder path for plots and dat files.
    %   data_filename   :   path to the .mat file containing raw simulation data.
    %   file_prefix     :   prefix to use for generated files.
    %
    %   Output:
    %   None. Generates .png and .dat files.
    %
    % ---------------------------------------------------------
    load(data_filename, 'resultsArray', 'experimentsParam', 'L_sim', 'L', 'delta', 'array_number_simulations');

    experimentsParam = build_experiment_functions(experimentsParam);

    target_pts = 1024;
    limit_plots = 8;
    limit_spectrogram = 4;

    omega_grid = resultsArray(1).omega_grid;
    S_phi_true = experimentsParam.M_k(0, omega_grid);
    [~, zero_idx] = min(abs(omega_grid));

    comparison_matrix_S = [omega_grid(:), S_phi_true(:)];

    % --- COMPUTE TRUE RHO_TF (rho_{1,Spec}) ---
    rhoTF_true = zeros(1, length(omega_grid));
    for i = 1:length(omega_grid)
        % The moments M_k are a numerical estimation of the integral $ t.^k .* exp(-2 * pi * (w - t).^2) .* S(t) dt$
        M0_val = experimentsParam.M_k(0, omega_grid(i));
        M1_val = experimentsParam.M_k(1, omega_grid(i));
        M2_val = experimentsParam.M_k(2, omega_grid(i));
        if M0_val < 1e-150
            rhoTF_true(i) = NaN;
        else
            % The formula below follows from expanding the derivatives in Eq. (4.1).
            rhoTF_true(i) = 4 * pi * (M2_val * M0_val - M1_val^2) / (M0_val^2);
        end
    end
    comparison_matrix_rho = [omega_grid(:), rhoTF_true(:)];

    for idx = 1:length(array_number_simulations)
        n_sims = array_number_simulations(idx);
        all_omega_zeros = resultsArray(idx).omega_zeros_all;
        all_tau_zeros = resultsArray(idx).tau_zeros_all;

        if ~isfield(resultsArray(idx), 'omega_zeros_by_realization') || isempty(resultsArray(idx).omega_zeros_by_realization)
            error('Missing per-realization data (omega_zeros_by_realization). Please re-run the simulations.');
        end

        omega_zeros_by_realization = resultsArray(idx).omega_zeros_by_realization;

        Z_all = zeros(n_sims, length(omega_grid));

        omega_pos = omega_grid(zero_idx:end);
        omega_neg = abs(omega_grid(zero_idx:-1:1));

        for k = 1:n_sims
            omega_zeros_k = omega_zeros_by_realization{k};

            % Zeros per frequency row of the grid.
            row = round((omega_zeros_k(:) - omega_grid(1)) / delta) + 1;
            row = row(row >= 1 & row <= length(omega_grid));
            n_row = accumarray(row, 1, [length(omega_grid) 1]).';

            density = n_row / delta;
            zero_count = zeros(1, length(omega_grid));
            zero_count(zero_idx:end)  = cumtrapz(omega_pos, density(zero_idx:end));
            zero_count(zero_idx:-1:1) = cumtrapz(omega_neg, density(zero_idx:-1:1));

            % Empirical integral Z_l(omega) of the zero-count.
            Z_l = zeros(1, length(omega_grid));

            Z_l(zero_idx:end) = cumtrapz(omega_pos, zero_count(zero_idx:end));
            Z_l(zero_idx:-1:1) = cumtrapz(omega_neg, zero_count(zero_idx:-1:1));

            Z_all(k, :) = (2 * pi / L) * Z_l;
        end

        Z_bar = mean(Z_all, 1);
        Z_var = var(Z_all, 0, 1);

        Z_var_of_mean = Z_var / n_sims;

        % Estimator of the smoothed PSD S_phi with log-normal bias correction (Eq. (5.29))
        S_phi_hat = exp(Z_bar - 2 * pi * (omega_grid.^2) - (Z_var_of_mean / 2));

        comparison_matrix_S = [comparison_matrix_S, S_phi_hat(:)];

        binSize = 0.2;
        offset = 0.005;
        edges = -limit_plots - offset:binSize:limit_plots - offset;
        centers = edges(1:end - 1) + binSize / 2;

        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
        if isOctave
            counts_raw = histc(all_omega_zeros, edges);
            counts = counts_raw(1:end - 1);
            counts(end) = counts(end) + counts_raw(end);
            counts = reshape(counts, 1, []);
        else
            counts = histcounts(all_omega_zeros, edges);
        end

        % Implements the estimator for the first intensity rho_{1,Spec} in (5.27)
        area_per_bin = binSize * 2 * L * n_sims;
        rhoTF_hat = counts / area_per_bin;

        % Mask overflowing values.
        clean_S_hat = S_phi_hat;
        clean_S_true = S_phi_true;
        clean_S_hat(~isfinite(clean_S_hat) | abs(clean_S_hat) > 1e30) = NaN;
        clean_S_true(~isfinite(clean_S_true) | abs(clean_S_true) > 1e30) = NaN;

        % Plot A: S_phi Estimation for a particular number of realizations
        figure(1200 + idx);
        clf;
        plot(omega_grid, clean_S_hat, 'r', 'LineWidth', 1.5);
        hold on;
        plot(omega_grid, clean_S_true, '--k', 'LineWidth', 2);
        hold off;
        grid on;
        xlim([-limit_plots, limit_plots]);

        valid_S = clean_S_true(~isnan(clean_S_true));
        if isempty(valid_S) || max(valid_S) <= 0
            max_S = 1;
        else
            max_S = max(valid_S);
        end
        ylim([0, max_S * 1.2]);

        xlabel('Frequency (\omega)', 'Interpreter', 'tex');
        hLeg = legend('Empirical S_{\phi}(\omega)', 'True S_{\phi}(\omega)');
        set(hLeg, 'Interpreter', 'tex');
        saveas(gcf, fullfile(folderPath, sprintf('%s_n%d_S_smoothed_individual.png', file_prefix, n_sims)));

        % Plot B: First intensity
        figure(1500 + idx);
        clf;
        bar(centers, rhoTF_hat, 1.0, 'FaceColor', [0.2, 0.2, 0.6], 'EdgeColor', 'none');
        hold on;
        plot(omega_grid, rhoTF_true, 'r', 'LineWidth', 2);
        hold off;

        xlabel('Frequency (\omega)', 'Interpreter', 'tex');
        hLeg = legend('Empirical \rho_{1,Spec}(\omega)', 'True \rho_{1,Spec}(\omega)');
        set(hLeg, 'Interpreter', 'tex', 'Location', 'northwest');
        xlim([-limit_plots, limit_plots]);

        valid_true = rhoTF_true(isfinite(rhoTF_true));
        valid_emp = rhoTF_hat(isfinite(rhoTF_hat));
        max_val = max([valid_true, valid_emp, 0.1]);
        ylim([0, max_val * 1.2]);

        saveas(gcf, fullfile(folderPath, sprintf('%s_n%d_rhoTF.png', file_prefix, n_sims)));

        % Plot C: First Spectrogram
        if idx == 1
            figure(1300 + idx);
            clf;
            global cRange
            cRange = [-100, -10];
            plotMatrix(abs(resultsArray(idx).first_stft.').^2, resultsArray(idx).tau_grid, resultsArray(idx).omega_grid);
            set(gca, 'YDir', 'normal');
            colorbar;
            hold on;
            scatter(resultsArray(idx).first_tau_zeros, resultsArray(idx).first_omega_zeros, 15, 'w', 'filled', 'MarkerEdgeColor', 'k');
            hold off;
            xlim([-limit_spectrogram, limit_spectrogram]);
            ylim([-limit_spectrogram, limit_spectrogram]);
            xlabel('Time (\tau)', 'Interpreter', 'tex');
            ylabel('Frequency (\omega)', 'Interpreter', 'tex');
            saveas(gcf, fullfile(folderPath, sprintf('%s_spectrogram.png', file_prefix)));
        end

        % Export individual .dat files
        disp(['Exporting individual .dat files for n = ', num2str(n_sims), '...']);
        n_total = length(omega_grid);
        if n_total > target_pts
            bin = min(floor((0:n_total - 1)' * target_pts / n_total) + 1, target_pts);
            cnt = accumarray(bin, 1, [target_pts 1]);

            export_omega   = accumarray(bin, omega_grid(:),      [target_pts 1]) ./ cnt;
            export_S_emp   = accumarray(bin, S_phi_hat(:),       [target_pts 1]) ./ cnt;
            export_S_true  = accumarray(bin, S_phi_true(:),      [target_pts 1]) ./ cnt;
            export_rhoTF_true = accumarray(bin, rhoTF_true(:),   [target_pts 1]) ./ cnt;
        else
            export_omega   = omega_grid(:);
            export_S_emp   = S_phi_hat(:);
            export_S_true  = S_phi_true(:);
            export_rhoTF_true = rhoTF_true(:);
        end
        mat_S = [export_omega, export_S_emp, export_S_true];
        save(fullfile(folderPath, sprintf('%s_n%d_S_smoothed.dat', file_prefix, n_sims)), 'mat_S', '-ascii');

        mat_rhoTF_emp = [centers(:), rhoTF_hat(:)];
        save(fullfile(folderPath, sprintf('%s_n%d_rhoTF_empirical.dat', file_prefix, n_sims)), 'mat_rhoTF_emp', '-ascii');

        mat_rhoTF_true_export = [export_omega, export_rhoTF_true];
        save(fullfile(folderPath, sprintf('%s_true_rhoTF_theoretical.dat', file_prefix)), 'mat_rhoTF_true_export', '-ascii');

    end

    % --- Combined Plot for S_phi for several realizations ---
    figure(1400);
    clf;
    plot(omega_grid, S_phi_true, '--k', 'LineWidth', 2.5);
    hold on;
    for i = 1:length(array_number_simulations)
        plot(omega_grid, comparison_matrix_S(:, i + 2), 'LineWidth', 1.2);
    end
    hold off;
    grid on;
    xlim([-limit_plots, limit_plots]);
    ylim([0, max(S_phi_true) * 1.2]);
    leg_str = cell(length(array_number_simulations) + 1, 1);
    leg_str{1} = 'True S_{\phi}';
    for i = 1:length(array_number_simulations)
        leg_str{i + 1} = sprintf('Estimator when n=%d', array_number_simulations(i));
    end
    hLeg = legend(leg_str);
    set(hLeg, 'Interpreter', 'tex');
    xlabel('Frequency (\omega)', 'Interpreter', 'tex');
    saveas(gcf, fullfile(folderPath, sprintf('%s_comparison_all_n_S_smoothed.png', file_prefix)));

    disp(['Exporting comparison .dat files for ', file_prefix, '...']);

    % Downsample S_phi comparison matrix for LaTeX
    n_total = size(comparison_matrix_S, 1);
    num_cols = size(comparison_matrix_S, 2);
    if n_total > target_pts
        bin = min(floor((0:n_total - 1)' * target_pts / n_total) + 1, target_pts);
        cnt = accumarray(bin, 1, [target_pts 1]);
        export_mat_S = zeros(target_pts, num_cols);
        for c = 1:num_cols
            export_mat_S(:, c) = accumarray(bin, comparison_matrix_S(:, c), [target_pts 1]) ./ cnt;
        end
    else
        export_mat_S = comparison_matrix_S;
    end
    save(fullfile(folderPath, sprintf('%s_comparison_S_smoothed.dat', file_prefix)), 'export_mat_S', '-ascii');

    disp('All .dat files generated successfully.');
end
