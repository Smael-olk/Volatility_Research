function MCSimulation()
%MCSIMULATION Verifies NIG Simulation algorithm across different maturities.
%   Loads parameters from 'CalibratedParams.mat' and compares MC prices
%   against a high-precision FFT benchmark for Short, Medium, and Long T.

    % --- 1. DATA LOADING & SETUP ---
    results_file = './params/CalibratedParams.mat';
    target_date  = '2017-12-08'; % Date to verify
    
    if ~isfile(results_file)
        error('File %s not found.', results_file);
    end
    fprintf('Loading calibration results...\n');
    load(results_file, 'Results');
    
    % Find the specific date index
    found_idx = -1;
    for i = 1:length(Results)
        if isequal(dateshift(Results(i).Date, 'start', 'day'), ...
                   dateshift(datetime(target_date), 'start', 'day'))
            found_idx = i;
            break;
        end
    end
    
    if found_idx == -1
        error('Date %s not found in results.', target_date);
    end
    
    % Extract Day Data
    DayData = Results(found_idx);
    
    % NIG Model Parameters (Calibrated)
    sigma_opt = DayData.Params(1);
    eta_opt   = DayData.Params(2);
    k_opt     = DayData.Params(3);
    
    fprintf('Analyzing Date: %s | Params: sigma=%.4f, eta=%.4f, k=%.4f\n', ...
            datestr(DayData.Date), sigma_opt, eta_opt, k_opt);

    % --- 2. CONFIGURATION FOR ANALYSIS ---
    % Select 3 distinct maturities to verify (Short, Middle, Long)
    % Using indices from the available term structure
    T_indices = [1, round(length(DayData.T_years)/2), length(DayData.T_years)];
    
    % Simulation Settings
    param = 0.01;       % FFT grid spacing (dz)
    N_sim = 10^6;       % Number of MC simulations (Reduced slightly for speed, incr. to 10^7 for final)
    alpha = 1/2;        % Standard NIG
    
    % Prepare Plot
    figure('Name', 'MC Convergence Across Maturities', 'Color', 'black','InvertHardcopy','off');
    hold on;
    colors = {'r', 'b', 'g'};
    markers = {'o', 's', '^'};
    styles = {'-', '--', '-.'};
    
    % --- 3. MAIN LOOP OVER MATURITIES ---
    for t_step = 1:length(T_indices)
        idx_T = T_indices(t_step);
        
        % Extract Market Variables for this Maturity
        Dt = DayData.T_years(idx_T);
        B  = DayData.Discounts(idx_T);
        F0 = DayData.Fwds(idx_T);
        
        fprintf('\n--- Testing Maturity T = %.2f years ---\n', Dt);
        
        % A. Define Characteristic Function (Specific to this Dt)
        % Calculate optimal shift 'a'
        p = (0.5 + eta_opt) + sqrt((0.5 + eta_opt)^2 + 2 * (1 - alpha) / (sigma_opt^2 * k_opt));
        a = p / 2; 
        
        % Define Laplace Exponent L(z)
        L = @(z, k, s) exp(Dt / k * (1 - sqrt(1 + 2 * k * z * s^2)));
        
        % Define CF phi(u)
        phi = @(u) exp(-1i * u * log(L(eta_opt, k_opt, sigma_opt))) .* ...
                   L(((u).^2 + 1i * u * (1 + 2 * eta_opt)) / 2, k_opt, sigma_opt);

        % B. Define Strikes (Around ATM)
        K = linspace(F0*0.9, F0*1.1, 40);
        
        % C. Calculate "TRUE" Benchmark Price (High Precision FFT)
        M_bench = 15; 
        BenchmarkPrices = PriceCallOption(K, F0, B, phi, M_bench, param, a);
        
        % D. Convergence Loop (Varying M from 4 to 14)
        M_values = 4:15;
        conv_error = zeros(size(M_values));
        
        for j = 1:length(M_values)
            m = M_values(j);
            
            % Simulate Samples
            samples = SimulateFromCF(phi, m, param, a, N_sim);
            
            % Price Option using MC
            % Payoff: max(S_T - K, 0)
            % S_T = F0 * exp(x_T)
            S_T = F0 * exp(samples);
            
            % Vectorized payoff calc for all strikes
            % Result is (1 x num_strikes)
            payoffs = max(S_T - K, 0); 
            MC_Price = mean(payoffs, 1) * B;
            
            % Compute Max Absolute Error vs Benchmark
            conv_error(j) = max(abs(MC_Price - BenchmarkPrices));
        end
        
        % E. Plot Convergence for this Maturity
        plot(M_values, log10(conv_error), ...
             'Color', colors{t_step}, ...
             'Marker', markers{t_step}, ...
             'LineStyle', styles{t_step}, ...
             'LineWidth', 2, ...
             'DisplayName', sprintf('T = %.2f yrs', Dt));
         
        fprintf('Max Error at M=14: %.2e\n', conv_error(end));
    end
    
    % --- 4. FORMAT PLOT ---
    xlabel('Grid Resolution M (N = 2^M)');
    ylabel('Log_{10}(Max Pricing Error)');
    title(['MC Algorithm Verification: ' datestr(DayData.Date)]);
    legend('Location', 'best');
    grid on;
    ylim([-5, 0]); % Adjust based on results
    hold off;
    
    disp('Verification Complete.');
end

% --- REQUIRED HELPER FUNCTIONS ---
% Ensure SimulateFromCF, PriceCallOption, etc., are in your path.