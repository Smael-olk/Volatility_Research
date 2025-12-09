function plot_nig_parameter_evo()
% PLOT_NIG_PARAMETER_EVOLUTION Loads calibrated NIG parameters and plots 
% their evolution over time to visualize stability and breakdown points.
    % --- 1. Load Data ---
    try
        % Load the structure containing the Results array
        NIG_data = load('./params/CalibratedParams.mat');
        Results = NIG_data.Results;
    catch
        error('Error loading CalibratedParams.mat. Ensure the file exists and contains the "Results" struct.');
    end
    
    if isempty(Results)
        disp('The Results structure is empty. Cannot plot.');
        return;
    end

    % --- 2. Extract Data ---
    
    % Initialize arrays for plotting
    numDays = length(Results);
    Dates = NaT(numDays, 1);
    Sigma = NaN(numDays, 1); % Volatility scaling parameter
    Eta = NaN(numDays, 1);   % Skewness parameter (asymmetry)
    Kappa = NaN(numDays, 1); % Kurtosis parameter (tail thickness, usually named 'k')
    
    % Loop through results to extract and clean data
    for j = 1:numDays
        % Check if parameters were successfully calibrated (not empty or NaN)
        if ~isempty(Results(j).Params) && ~any(isnan(Results(j).Params))
            Dates(j) = Results(j).Date;
            Sigma(j) = Results(j).Params(1);
            Eta(j) = Results(j).Params(2);
            Kappa(j) = Results(j).Params(3);
        end
    end

    % --- 3. Clean and Filter Data ---
    
    % Remove days where calibration failed (Date is NaT or parameters are NaN)
    valid_indices = ~isnat(Dates) & ~isnan(Sigma);
    Dates = Dates(valid_indices);
    Sigma = Sigma(valid_indices);
    Eta = Eta(valid_indices);
    Kappa = Kappa(valid_indices);

    if isempty(Dates)
        disp('No valid calibrated data points found.');
        return;
    end

    % --- 4. Plotting ---
    
    figure('Name', 'NIG Parameter Evolution', 'Position', [100 100 1200 800]);
    h_fig = gcf;
    set(h_fig, 'InvertHardcopy', 'off');
    set(h_fig, 'Renderer', 'painters');

    % Subplot 1: Sigma (Volatility Scaling)
    subplot(3, 1, 1);
    plot(Dates, Sigma, 'b-', 'LineWidth', 2);
    title('1. Volatility Scaling Parameter (\sigma)');
    ylabel('\sigma Value');
    grid on;

    % Subplot 2: Eta (Skewness / Asymmetry)
    subplot(3, 1, 2);
    plot(Dates, Eta, 'r-', 'LineWidth', 2);
    title('2. Skewness Parameter (\eta)');
    ylabel('\eta Value');
    grid on;

    % Subplot 3: Kappa (Kurtosis / Tail Thickness)
    subplot(3, 1, 3);
    plot(Dates, Kappa, 'g-', 'LineWidth', 2);
    title('3. Kurtosis Parameter (\kappa)');
    ylabel('\kappa Value');
    xlabel('Trade Date');
    grid on;
    
    % Add a super title for clarity
    sgtitle('Time Series Evolution of Calibrated NIG Parameters', 'FontWeight', 'bold');

    % --- 5. Identify Breakdown Period ---
    % A breakdown often appears as a sharp spike, a large jump, or a flat line 
    % where the optimizer hit a boundary. Use the plots to find the date range.
    fprintf('\nAnalysis Complete. Check plots for:\n');
    fprintf(' - Sharp Spikes/Jumps: Indicate parameter instability.\n');
    fprintf(' - Flat Lines: Suggest optimizer hitting a bound (e.g., lower bound of sigma/kappa).\n');
    fprintf(' - Erratic Movements: Sign of local minima traps/non-convergence.\n');
end