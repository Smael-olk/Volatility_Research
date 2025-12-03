function mkt_vol_surface_evo()
%VOL_SURFACE_EVO Summary of this function goes here
%   Detailed explanation goes here

    filelist = dir('.\Dati_Train\*.csv'); 
    numFiles = length(filelist);
    spot_file_name = 'spot_SPX_hist.csv'; % Name of the spot price file
    main_fig = figure('Name', 'Dynamic Volatility Surface', 'Position', [100 100 1200 600]);
    
    % Set default figure properties for black background and white text
    set(0, 'DefaultFigureColor', [0 0 0]); % Black background
    set(0, 'DefaultAxesColor', [0 0 0]);   % Black axes background
    set(0, 'DefaultAxesXColor', [1 1 1]);  % White X-axis labels/ticks
    set(0, 'DefaultAxesYColor', [1 1 1]);  % White Y-axis labels/ticks
    set(0, 'DefaultTextColor', [1 1 1]);   % white text 
    main_ax = axes('Parent', main_fig);
    
    v = VideoWriter('mkt_vol_surface_Animation.avi');
    v.FrameRate = 7; 
    open(v);
    for j = 1:numFiles
        current_options_file = filelist(j).name;
        [discounts,forward_prices,expiries,TradeDate,Spot,Option_data]=IR_curves(current_options_file,spot_file_name); % Discounts and Forward Calibration function
        date_cell_array = cellstr(datestr(expiries));
        dates_str = strjoin(date_cell_array, ', ');
        T_years = years(expiries - TradeDate);
        
        % --- STEP C: Plotting and Animation Capture ---
        if ~isempty(T_years)
            marketdata = cell(length(Option_data.datesExpiry), 8);
            for i=1:length(Option_data.datesExpiry)
                marketdata{i,1} = T_years(i);
                marketdata{i,2} = discounts(i);
                marketdata{i,3} = forward_prices(i);
                marketdata{i,4} = (Option_data.callAsk(i).prices + Option_data.callBid(i).prices ) / 2; % observed market prices
                marketdata{i,5} = Option_data.strikes(i).value;
            end
            
            sigma = 0.1; % Initial guess (often 0.2 to 0.4 works well)
            tol = 1e-6;
            max_iter = 200;
            for i=1:length(Option_data.datesExpiry)
                iv= [];
                Dt = marketdata{i,1};
                B= marketdata{i,2};
                F= marketdata{i,3};
                K_market = marketdata{i,5};
                C=marketdata{i,4};
                for j=1:length(C)
                    iv(j)= Black76Inverse(C(j), F, K_market(j), B, Dt, tol, max_iter);
                end
                marketdata{i,7} = iv;
            end
            
            vol_surface_plot(marketdata, main_ax,TradeDate);

            % Capture the frame for the video
            frame = getframe(main_fig);
            writeVideo(v, frame);
        end
    end
    
    close(v);
    disp('Animation saved as vol_surface_Animation.avi');
end

