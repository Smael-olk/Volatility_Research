function  term_structure_evo()
%TERM_STRUCTURE_EVO Summary of this function goes here
%   Detailed explanation goes here

    filelist = dir('.\Dati_Train\*.csv'); 
    numFiles = length(filelist);
    spot_file_name = 'spot_SPX_hist.csv'; % Name of the spot price file
    figure('Name', 'Dynamic Term Structure', 'Position', [100 100 1200 600]);
    
    % Set default figure properties for black background and white text
    set(0, 'DefaultFigureColor', [0 0 0]); % Black background
    set(0, 'DefaultAxesColor', [0 0 0]);   % Black axes background
    set(0, 'DefaultAxesXColor', [1 1 1]);  % White X-axis labels/ticks
    set(0, 'DefaultAxesYColor', [1 1 1]);  % White Y-axis labels/ticks
    set(0, 'DefaultTextColor', [1 1 1]);   % white text 
    
    
    figure('Name', 'Dynamic Term Structure', 'Position', [100 100 1200 600]);
    v = VideoWriter('./Plots/Term_Structure_Animation.avi');
    v.FrameRate = 2; 
    open(v);
    for j = 1:numFiles
        current_options_file = filelist(j).name;
        [discounts,forward_prices,expiries,TradeDate,Spot]=IR_curves(current_options_file,spot_file_name,0.1, 0.6, 5, 50); % Discounts and Forward Calibration function
        date_cell_array = cellstr(datestr(expiries));
        dates_str = strjoin(date_cell_array, ', ');
        T_years = years(expiries - TradeDate);
        
        % --- STEP C: Plotting and Animation Capture ---
        if ~isempty(T_years)
            clf; % Clear the figure for the new plot
    
            % Define the custom orange color (RGB values)
            orange_color = [251/255 139/255 30/255];
    
            % 1. Zero Rate Curve (Left Subplot)
            R_years = -log(discounts(:)') ./ T_years(:)'; % Calculate the Zero Rate
            
            subplot(1, 2, 1);
            plot(T_years, R_years * 100, 'Color', orange_color, 'LineWidth', 2, 'MarkerSize', 6, 'Marker', '^');
            
            title(['Zero Rate Curve: ', datestr(TradeDate, 'yyyy-mm-dd')]);
            xlabel('Time to Maturity (Years)');
            ylabel('Zero Rate r(T) (%)');
            
            % Ensure grid lines are visible on black background (optional, but good practice)
            set(gca, 'GridColor', [0.3 0.3 0.3]); % Dark gray grid for visibility
            grid on;
            axis([0 max(T_years)+0.5 0.5 3]); 
            ylim([0,5])
        
            % 2. Forward Price Curve (Right Subplot)
            subplot(1, 2, 2);
            plot(T_years, forward_prices, 'Color', orange_color, 'LineWidth', 2, 'MarkerSize', 6, 'Marker', 's');
            
            title(['Forward Price Curve: Spot=', num2str(Spot, '%.2f')]);
            xlabel('Time to Maturity (Years)');
            ylabel('Forward Price F(t_0, T)');
            
            set(gca, 'GridColor', [0.3 0.3 0.3]); % Dark gray grid
            grid on;
            
            axis([0 max(T_years)+0.5 Spot-10 Spot+100]); 
            ylim([Spot*0.96,Spot*1.05])
            % Capture the frame for the video
            frame = getframe(gcf);
            writeVideo(v, frame);
        end
    end
    
    close(v);
    disp('Animation saved as Term_Structure_Animation.avi');
end

