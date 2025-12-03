%% 1. Initialization and Configuration
clear; close all; clc;

% --- Simulation Settings ---
flag_SP = 1;                % 1 for S&P 500, 0 for EUROSTOXX
DK = 5;                     % Minimum difference between strikes in the dataset

% --- Date Settings ---
start_date = datetime('2017-06-08', 'InputFormat', 'yyyy-MM-dd');
stop_date  = datetime('2017-12-08', 'InputFormat', 'yyyy-MM-dd');

target_maturity = datetime('19-Jan-2018'); % Maturity of the option we are hedging

% --- Physics & Thresholds ---
weekend_mask = [1, 0, 0, 0, 0, 0, 1]; % 1 = Weekend
delta = 0.1;        % Delta Trheshold        
tresh_bid_ask = 0.6;
tresh_penny = 0.1;
thresh_volume = 0; %keeping al value
x0_sigma = 0.2;             % Initial sigma guess (For BLS inv)


% --- FFT Param ---
calib_param = 0.01;
M = 14;                     % FFT sample parameter

% --- File Paths ---
% Use fullfile for robust path handling (avoids issues with '\' vs '/')
base_location = '.\';

if flag_SP
    % S&P 500 Configuration
    surf_name = "S&P";
    sub_dir = fullfile('Dati_Train');
    
    namefile_spot = fullfile(base_location, 'spot_SPX_hist.csv');
    
else
    % EUROSTOXX Configuration
    surf_name = "EU";
    sub_dir = fullfile('SurfaceStorica', 'Surface_Giornaliera_EUR');
    
    namefile_spot = fullfile(base_location, sub_dir, 'spot_EUR_hist.csv');
    
end

full_sub_path = fullfile(base_location, sub_dir);

%% 2. Main Loop Setup
% Pre-allocation or initialization of history vectors
price_vec = [];
price_hedged_vec = [];
F_hist = [];
ix_act = [];
NIG_mat=[];
date_vec=[];
% Portfolio State Initialization
Spot_to_hedge_old = 0;
cash_old = 0;
iter = start_date;

disp(['Starting simulation for: ', char(surf_name)]);

%% 3. Simulation Loop
date_range = start_date:stop_date;
r_old=0.04;
refDate_old=start_date;
refDate = datetime('08-Dec-2019');

for i = 1:length(date_range)
    current_date = date_range(i);
    
    % Construct daily filename
    date_str = datestr(current_date, 'yyyy-mm-dd');
    filename_opt = [date_str, '.csv'];
    path_opt = fullfile(full_sub_path, filename_opt);
    
    % Check if Business Day AND File Exists
    is_bus_day = (business_day(current_date, weekend_mask) == current_date);
    
    if is_bus_day && isfile(path_opt)
        date_vec=[date_vec,current_date];
        % --- A. Data Loading ---
        %namefile_lib = fullfile(full_sub_path, 'LIB_hist.csv');
        
           
        % Option Extraction
        mkt = Option_estraction_hist(path_opt, namefile_spot, tresh_bid_ask, tresh_penny, DK);
        Ft = mkt.spot*exp(-yearfrac(refDate_old,refDate)*r_old);
    
        %Calibrate
        % --- C. Portfolio Evaluation ---
        % Time to maturity
        Dt = yearfrac(current_date, target_maturity, 3); 
        
      
        
        % Calculate Derivative Price (Long call option)
        %[Price_port] = evaluate_portfolio(Param_NIG, B_vec, Ft, maturity, Dt, calib_param, M);
        for i=1:length(mkt.strikes)
            [Price_port] = blsprice(mkt.spot,mkt.strikes(i).value,r_old,Dt,x0_sigma);
            %Calculate IR remuneration for the period
            B_old=exp(-yearfrac(refDate_old,refDate)*r_old);
            refDate_old=refDate;
            % Value of the hedged portfolio
            Port_value_hedge = (Spot_to_hedge_old * mkt.spot) + Price_port + cash_old/B_old;
        end
        % --- D. Hedging Strategy (Delta Hedging) ---
        % Finite difference calculation for Delta
        %price_bumped = evaluate_portfolio(Param_NIG, B_vec, Ft + 1, maturity, Dt, calib_param, M);
        price_bumped =  blsprice(mkt.spot+1,mkt.strikes(i).value,r_old,Dt,x0_sigma);
        delta_proxy = (price_bumped - Price_port) * (Ft / mkt.spot); 
        Spot_to_hedge = -delta_proxy;
        
        % Update Cash Account
        cash = (Spot_to_hedge_old - Spot_to_hedge) * mkt.spot + cash_old/B_old;
        
        % Update State Variables
        Spot_to_hedge_old = Spot_to_hedge;
        cash_old = cash;
        Parameters_old = Parameters;
        r_old=-log(B(1))/yearfrac(refDate,dates(1),3);
     
        % Store History 
        price_vec = [price_vec, Price_port];
        price_hedged_vec = [price_hedged_vec, Port_value_hedge]; 
        F_hist = [F_hist, Ft]; 
        NIG_mat=[NIG_mat,Param_NIG]
        % Progress indicator (optional)
        fprintf('Processed: %s\n', date_str);

    end
end

%% 4. Visualization
figure('Name', 'Hedge Simulation Results', 'NumberTitle', 'off');
hold on; grid on;

% Normalize to zero at start for comparison
plot(price_vec - price_vec(1), 'LineWidth', 1.5);
plot(F_hist - F_hist(1), 'LineWidth', 1.5);
plot(price_hedged_vec - price_hedged_vec(1), 'LineWidth', 1.5, 'LineStyle', '--');

title(['Hedging Performance: ', char(surf_name)]);
xlabel('Trading Days');
ylabel('Relative Value');
legend("Derivative Price","Forward Price","Hedged Portfolio", 'Location', 'best');

% Normalize to zero at start for comparison
figure
hold on
plot( diff(price_vec).^2, 'LineWidth', 1.5);
plot(diff(price_hedged_vec).^2, 'LineWidth', 1.5, 'LineStyle', '--');

title(['Hedging Performance (Squared P&L) ', char(surf_name)]);
xlabel('Trading Days');
ylabel('Relative Value');
legend("Derivative Price", "Hedged Portfolio", 'Location', 'best');

hold off;

figure
title(['Evolution of sigma', char(surf_name)],'LineWidth', 1.5);

plot(NIG_mat(2:3:end))

figure
title(['Evolution of eta', char(surf_name)],'LineWidth', 1.5);

plot(NIG_mat(3:3:end))


figure
title(['Evolution of k', char(surf_name)],'LineWidth', 1.5);

plot(NIG_mat(1:3:end))