function  Hedging(StartDate,EndDate,FinalDate)
%HEDGING Summary of this function goes here

Results_File = './params/CalibratedParams.mat';
% Positions (Initial)
Pos_Cert = -1;       % Short the Certificate
Pos_Underlying = 0; 
Pos_Put = 0;         % Risk Reversal Leg 1 (Short Put)
Pos_Call = 0;        % Risk Reversal Leg 2 (Long Call)
Cash_Account = 0; 

% Storage
Prev_Port_Value = 0;
PnL_History = [];
Dates_History = [];

% --- 2. LOAD & FILTER DATA ---
if ~isfile(Results_File), error('File not found: %s', Results_File); end
load(Results_File, 'Results');

% Extract dates and sort to ensure chronological order
[AllDates, sortIdx] = sort([Results.Date]);
SortedResults = Results(sortIdx);

% Filter for the 6-month window
InRange = (AllDates >= StartDate) & (AllDates <= EndDate);
BacktestData = SortedResults(InRange);

if isempty(BacktestData)
    error('No data found between %s and %s', StartDate, EndDate);
end

fprintf('Starting Backtest: SHORT RISK REVERSAL (Delta-Vega-Gamma)\n');
fprintf('Period: %s to %s (%d Trading Days)\n', ...
    datestr(StartDate), datestr(EndDate), length(BacktestData));

% --- 3. MAIN LOOP ---
for t = 1:length(BacktestData)
    DayData = BacktestData(t);
    CurrentDate = DayData.Date;
    
    params  = DayData.Params; % [alpha, beta, delta, mu, etc.]
    S0      = DayData.Spot;
    
    % Interpolate Rates (Yield Curve)
    T_final = years(FinalDate - CurrentDate);
    r_curve = -log(DayData.Discounts) ./ DayData.T_years;
    
    % Handle rate extrapolation safely
    r_riskfree = interp1(DayData.T_years, r_curve, T_final, 'linear', 'extrap');
    F_2y = interp1(DayData.T_years, DayData.Fwds, T_final, 'linear', 'extrap');
    q_div = r_riskfree - log(F_2y / S0) / T_final;
    
    % --- 4. VALUATION (Mark-to-Market) ---
    % Simulating existence of the Dec 2017 product back in June 2017
    [Price_Cert, ~] = Price_Miami_Certificate_Function(params, S0, r_riskfree, q_div, CurrentDate);
    Val_Cert = Pos_Cert * Price_Cert;
        
    % Value Existing Options (if any)
    Val_Options = 0;
    
    % -- Value Put --
    if Pos_Put ~= 0
        T_rem = years(Current_T - CurrentDate);
        if T_rem > 0.005 % If not expired
            P_Put_Curr = NIG_Pricer(params, S0, K_Put_Curr, exp(-r_riskfree*T_rem), T_rem, 2); 
            Val_Options = Val_Options + (Pos_Put * P_Put_Curr);
        else % Expired
            Payoff = max(K_Put_Curr - S0, 0);
            Cash_Account = Cash_Account + (Pos_Put * Payoff);
            Pos_Put = 0;
        end
    end
    
    % -- Value Call --
    if Pos_Call ~= 0
        T_rem = years(Current_T - CurrentDate);
        if T_rem > 0.005
            P_Call_Curr = NIG_Pricer(params, S0, K_Call_Curr, exp(-r_riskfree*T_rem), T_rem, 1);
            Val_Options = Val_Options + (Pos_Call * P_Call_Curr);
        else
            Payoff = max(S0 - K_Call_Curr, 0);
            Cash_Account = Cash_Account + (Pos_Call * Payoff);
            Pos_Call = 0;
        end
    end

    Val_Stock = Pos_Underlying * S0;
    
    % --- 5. CALCULATE P&L ---
    if t > 1
        % Interest on Cash
        dt_days = days(CurrentDate - BacktestData(t-1).Date);
        Interest = Cash_Account * (r_riskfree * dt_days/365);
        Cash_Account = Cash_Account + Interest;
        
        Port_Value = Val_Cert + Val_Stock + Val_Options + Cash_Account;
        Daily_PnL = Port_Value - Prev_Port_Value;
        PnL_History(end+1) = Daily_PnL;
        Dates_History = [Dates_History; CurrentDate];
        
        % Optional: Less verbose printing for long loop
        if mod(t, 20) == 0 || t == length(BacktestData)
             fprintf('Date: %s | P&L: $%7.2f | Spot: %.2f\n', datestr(CurrentDate), Daily_PnL, S0);
        end
    else
        % Initialization Day (June 8)
        Cash_Account = -Val_Cert; % Fund the initial short
        Prev_Port_Value = 0; 
        fprintf('Initialized Portfolio on %s at Spot %.2f\n', datestr(CurrentDate), S0);
    end
    Prev_Port_Value = Val_Cert + Val_Stock + Val_Options + Cash_Account;
    
    % --- 6. CALCULATE GREEKS ---
    dS = S0 * 0.01; dVol = params(1) * 0.01;
    Px = @(s, p) Price_Miami_Certificate_Function(p, s, r_riskfree, q_div, CurrentDate);
    
    % Cert Greeks
    P_0 = Px(S0, params);
    Delta_Cert = (Px(S0+dS, params) - Px(S0-dS, params)) / (2*dS);
    Gamma_Cert = (Px(S0+dS, params) - 2*P_0 + Px(S0-dS, params)) / (dS^2);
    
    p_up = params; p_up(1) = params(1) + dVol;
    p_dn = params; p_dn(1) = params(1) - dVol;
    Vega_Cert  = (Px(S0, p_up) - Px(S0, p_dn)) / (2*dVol);
    
    Risk_Vega  = Pos_Cert * Vega_Cert;
    Risk_Gamma = Pos_Cert * Gamma_Cert;
    Risk_Delta = Pos_Cert * Delta_Cert;
    
    % --- 7. RE-HEDGE (Risk Reversal) ---
    % Strategy: 3-Month Options (Constant Maturity Roll)
    T_new = 0.25; 
    Target_Exp = CurrentDate + days(90);
    
    % Strikes: 95% Put, 105% Call (Low Delta / OTM)
    K_New_Put  = S0 * 0.95;
    K_New_Call = S0 * 1.05;
    
    % Helper: 1=Call, 2=Put
    OptPx = @(s, p, k, type) NIG_Pricer(p, s, k, exp(-r_riskfree*T_new), T_new, type);
    
    % Put Greeks
    P_Put0 = OptPx(S0, params, K_New_Put, 2);
    D_Put  = (OptPx(S0+dS, params, K_New_Put, 2) - OptPx(S0-dS, params, K_New_Put, 2)) / (2*dS);
    G_Put  = (OptPx(S0+dS, params, K_New_Put, 2) - 2*P_Put0 + OptPx(S0-dS, params, K_New_Put, 2)) / (dS^2);
    V_Put  = (OptPx(S0, p_up, K_New_Put, 2) - OptPx(S0, p_dn, K_New_Put, 2)) / (2*dVol);
    
    % Call Greeks
    P_Call0 = OptPx(S0, params, K_New_Call, 1);
    D_Call  = (OptPx(S0+dS, params, K_New_Call, 1) - OptPx(S0-dS, params, K_New_Call, 1)) / (2*dS);
    G_Call  = (OptPx(S0+dS, params, K_New_Call, 1) - 2*P_Call0 + OptPx(S0-dS, params, K_New_Call, 1)) / (dS^2);
    V_Call  = (OptPx(S0, p_up, K_New_Call, 1) - OptPx(S0, p_dn, K_New_Call, 1)) / (2*dVol);

    % Solve for Weights (Kill Vega & Gamma)
    A = [V_Put, V_Call; G_Put, G_Call];
    B = [-Risk_Vega; -Risk_Gamma];
    
    % Protect against singular matrix (rare, but good practice)
    if rcond(A) < 1e-10
        warning('Matrix singular at %s. Keeping previous hedge.', datestr(CurrentDate));
        Target_Pos_Put = Pos_Put; Target_Pos_Call = Pos_Call; Target_Pos_Stock = Pos_Underlying;
    else
        W_Opts = A \ B; 
        Target_Pos_Put  = W_Opts(1);
        Target_Pos_Call = W_Opts(2);
        
        Total_Opt_Delta = Target_Pos_Put * D_Put + Target_Pos_Call * D_Call;
        Target_Pos_Stock = -(Risk_Delta + Total_Opt_Delta);
    end

    % --- 8. EXECUTE TRADES ---
    % Close Old
    if Pos_Put ~= 0,  Cash_Account = Cash_Account - (-Pos_Put * P_Put_Curr); end
    if Pos_Call ~= 0, Cash_Account = Cash_Account - (-Pos_Call * P_Call_Curr); end
    
    % Open New
    Cost_New_Put  = Target_Pos_Put  * P_Put0;
    Cost_New_Call = Target_Pos_Call * P_Call0;
    Cash_Account = Cash_Account - Cost_New_Put - Cost_New_Call;
    
    % Adjust Stock
    Trade_Stock = Target_Pos_Stock - Pos_Underlying;
    Cash_Account = Cash_Account - (Trade_Stock * S0);
    
    % Update Portfolio
    Pos_Put = Target_Pos_Put;
    Pos_Call = Target_Pos_Call;
    Pos_Underlying = Target_Pos_Stock;
    K_Put_Curr = K_New_Put;
    K_Call_Curr = K_New_Call;
    Current_T = Target_Exp;
end

% Summary Stats
MSE = mean(PnL_History.^2);
Cum_PnL = sum(PnL_History);
fprintf('\n--------------------------------------------------\n');
fprintf('BACKTEST COMPLETE\n');
fprintf('Cumulative P&L: $%.2f\n', Cum_PnL);
fprintf('Mean Squared Error (MSE): %.4f\n', MSE);
% Optional: Plot P&L
bbg_orange = [1.0 0.6 0.0];
figure; plot(Dates_History, PnL_History,'color',bbg_orange); 
title('Daily Hedging P&L (June - Dec)'); xlabel('Date'); ylabel('P&L ($)');
grid on;
end

