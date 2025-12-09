function vol_surface_evo()
%VOL_SURFACE_EVO Generates video animation of the volatility surface 
%                from pre-calculated calibration results.

    % --- 1. Load Data ---
    if ~isfile('./params/CalibratedParams.mat')
        error('Results file not found. Run run_calibration.m first.');
    end
    load('./params/CalibratedParams.mat', 'Results'); 
    
    numFiles = length(Results);
    
    % --- 2. Figure and Video Setup ---
    main_fig = figure('Name', 'Dynamic Volatility Surface', 'Position', [100 100 1200 600]);
    
    % Set default figure properties for visualization
    set(0, 'DefaultFigureColor', [0 0 0]); 
    set(0, 'DefaultAxesColor', [0 0 0]);
    set(0, 'DefaultAxesXColor', [1 1 1]);
    set(0, 'DefaultAxesYColor', [1 1 1]);
    set(0, 'DefaultTextColor', [1 1 1]); 
    
    main_ax = axes('Parent', main_fig); % Axes handle for plotting
    
    % --- HIGH QUALITY VIDEO WRITER SETUP ---
    v = VideoWriter('./Plots/vol_surface_Animation_HQ.avi', 'Motion JPEG AVI'); 
    v.Quality = 99; % High quality
    v.FrameRate = 2; % 4 frames per second
    open(v);
    
    disp('Starting Animation Generation...');
    
    % --- 3. Plotting Loop ---
    % Assume parameters for Black76Inverse are stable:
    tol = 1e-6; 
    max_iter = 200;
    
    for j = 1:numFiles
        
        % Skip if calibration failed for this day (params is empty)
        if isempty(Results(j).Params)
            continue; 
        end
        
        % Unpack data from struct
        optimal_params = Results(j).Params;
        Spot = Results(j).Spot;
        TradeDate = Results(j).Date;
        T_years = Results(j).T_years;
        discounts = Results(j).Discounts;
        fwds = Results(j).Fwds;
        
        plot_data = cell(length(T_years), 8);
        
        % --- REGENERATE DENSE GRID PRICES AND IVs ---
        % This ensures a smooth IV surface for plotting, using the calibrated model
        for i = 1:length(T_years)
            Dt = T_years(i);
            B = discounts(i);
            F = fwds(i);
            
            plot_data{i,1} = Dt; % T
            
            % Create Dense Strike Grid (30 points for plotting smoothness)
            K_dense = linspace(Spot*0.85, Spot*1.2, 50); 
            plot_data{i,8} = K_dense; % K_dense
            
            % Calculate Model Prices (Flag 1 for Call)
            C_model = NIG_Pricer(optimal_params, F, K_dense, B, Dt);
            
            % Invert Model Prices to get Model-Consistent IVs
            iv = zeros(size(C_model));
            for k = 1:length(C_model)
                iv(k) = Black76Inverse(C_model(k), F, K_dense(k), B, Dt, tol, max_iter);
            end
            plot_data{i,7} = iv; % IVs
        end
        
        % --- 4. Plotting and Capture ---
        vol_surface_plot(plot_data, main_ax, TradeDate, Spot);
        
        % Capture the frame using the figure handle
        frame = getframe(main_fig); 
        writeVideo(v, frame);
    end
    
    close(v);
    disp('Animation saved as vol_surface_Animation_HQ.avi');
end