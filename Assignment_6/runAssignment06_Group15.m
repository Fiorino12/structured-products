% runAssignment06_Group15
%  Group 15, AY2023-2024
%
% to run:
% > runAssignment06_Group15
% Ferrari Irene, Fioravanti Mattia, Serva Giorgio Enrico Maria

clear 
close all;  
clc;

%% Settings and Loadings

% Data settings

formatData='dd/mm/yyyy';
zRatesconvention = 3; % Act/365
capletconvention = 2; % Act/360

%% Read market data

[datesSet, ratesSet] = readExcelData_bootstrap('MktData_CurveBootstrap_20-2-24', formatData);
Data_capvol = readExcelData_capvol( 'Caps_vol_20-2-24', datesSet.settlement);

%% Bootstrap discount factors

% Do the bootstrap of discount factors
[dates, discounts] = bootstrap(datesSet, ratesSet);
dates_yfz = yearfrac(dates(1), dates,zRatesconvention);
dates_yfc = yearfrac(dates(1), dates,capletconvention);
zRates = zeroRates(dates,discounts)/100;

% Plot discount and zero-rates curves
figure
date_plot = datetime(dates, 'ConvertFrom', 'datenum','Format','dd/MM/uuuu');
yyaxis left
plot(date_plot, discounts, '-s', 'MarkerFaceColor', 'green', 'MarkerSize', 5);
ylabel('discount factors');
xlabel('dates');
ylim([0.30,1])
hold on 
yyaxis right
ylabel('zero rates');
plot(date_plot, zRates,'-s', 'MarkerFaceColor', 'blue', 'MarkerSize', 5);
xlim([date_plot(1),date_plot(end)])
title('Bootstrap Curve')
grid
legend('discounts', 'Zero rates')
set(gca, 'XTick', date_plot , 'XTickLabelRotation',90)
set(gca,'XGrid','off','YGrid','on')
xticks([date_plot(1:6:end)])

%% Bootstrap volatilities

% Compute Cap market prices via LMM model with flat volatilities
Data_capvol.cap_prices = Price_Cap_flat(Data_capvol, dates, zRates); 

% Compute spot volatilities
Data_capvol.sigma_spot = bootstap_vol(Data_capvol,dates, zRates);

% Plot the flat and spot volatility surfaces
date_expiry_yfc = yearfrac(dates(1), finddates(dates(1), Data_capvol.expyear), capletconvention);
[Strike_spot_vol,Expiry_spot_vol] = meshgrid(Data_capvol.strikes,date_expiry_yfc);
% First subplot for Spot Volatilities
figure
subplot(1, 2, 1);
surf(Strike_spot_vol, Expiry_spot_vol, Data_capvol.sigma_spot);
title('Surface Plot of Spot Volatilities');
xlabel('Strike');
ylabel('Expiry');
zlabel('Spot Volatilities');
colorbar;
shading interp;
view(-30, 30);
grid on;
colormap(jet); 

% Second subplot for Flat Volatilities
subplot(1, 2, 2);
surf(Strike_spot_vol, Expiry_spot_vol, Data_capvol.flat_volatilities);
title('Surface Plot of Flat Volatilities');
xlabel('Strike');
ylabel('Expiry');
zlabel('Flat Volatilities');
colorbar;
shading interp;
view(-60, 40);
grid on;
colormap(jet); 


%% Pricing Certificate

% Import the certificate data 
Certificate_data = struct(...
                    "Notional" , 50*1e6,...
                    "T", 15,...
                    "spol_A", 2/100,...
                    "spol_B" , 1.1/100,...
                    "first_libor" , 3/100,...
                    "strike_5" , 4.3/100,...
                    "strike_5_10" , 4.6/100,...
                    "strike_10_15" , 5.1/100);


fprintf("------------------- Certificate upfront -------------------\n\n");

% Compute the certificate upfront
X = certificate_upfront(Certificate_data, dates, zRates,Data_capvol);

fprintf(" The certificate has an upfront of %.2f %% \n\n", X*100);


%% Compute the bucket sensitivities

% Define the elements used in the bootstrap 
n_depos = 4;
n_futures = 7;
n_swaps = 18;

% Initialize the vector of the bucket sensitivities
bucket_DV01 = zeros(n_depos+n_futures+n_swaps,1);

% Swap useful for hedging in next points
swap_dates = finddates(dates(1), (1:15)');
swap_rates =  mean(ratesSet.swaps([2,5,10,13],:),2);
bucket_swap_DV01 = zeros(n_depos+n_futures+n_swaps,4);

fprintf("------------------------------------------------------------ \n\n");

fprintf(" Starting computation of delta bucket sensitivities \n");

% Set the increment of 
h = 1e-8;
tic

% Compute the bucket delta sensitivites for all depos
for ii = 1:n_depos
    
    % Increase the single mkt rate
    ratesSet_DV01 = ratesSet;
    ratesSet_DV01.depos(ii,:) = ratesSet_DV01.depos(ii,:) + h;

    % Do the bootstrap with new mkt rates
    [~, discounts_DV01] = bootstrap(datesSet, ratesSet_DV01);
    zRates_DV01 = zeroRates(dates,discounts_DV01)/100;

    % Compute the new Cap market prices and recalibrate the spot vol surface 
    Data_capvol_DV01 = Data_capvol;
    Data_capvol_DV01.cap_prices = Price_Cap_flat(Data_capvol_DV01, dates, zRates_DV01); 
    Data_capvol_DV01.sigma_spot = bootstap_vol(Data_capvol_DV01,dates, zRates_DV01);
    
    % Price the certificate with the shifted data 
    X_bucket = certificate_upfront(Certificate_data, dates, zRates_DV01,Data_capvol_DV01);
    
    % Compute the certificate bucket sensitivity 
    bucket_DV01(ii) = (X_bucket - X)/h * 1e-4; 

    % Compute the swap bucket sensitivity 
    [bucket_swap_DV01(ii,1),~, ~] = sensSwap(dates(1), swap_dates(1:2), swap_rates(1), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(ii,2),~, ~] = sensSwap(dates(1), swap_dates(1:5), swap_rates(2), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(ii,3),~, ~] = sensSwap(dates(1), swap_dates(1:10), swap_rates(3), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(ii,4),~, ~] = sensSwap(dates(1), swap_dates(1:15), swap_rates(4), dates, discounts,discounts_DV01);
end


% Compute the bucket delta sensitivites for all futures
for ii = 1:n_futures
    % Increase the single mkt rate
    ratesSet_DV01 = ratesSet;
    ratesSet_DV01.futures(ii,:) = ratesSet_DV01.futures(ii,:) + h;

    % Do the bootstrap with new mkt rates
    [~, discounts_DV01] = bootstrap(datesSet, ratesSet_DV01);
    zRates_DV01 = zeroRates(dates,discounts_DV01)/100;
    
    % Compute the new Cap market prices and recalibrate the spot vol surface 
    Data_capvol_DV01 = Data_capvol;
    Data_capvol_DV01.cap_prices = Price_Cap_flat(Data_capvol_DV01, dates, zRates_DV01); 
    Data_capvol_DV01.sigma_spot = bootstap_vol(Data_capvol_DV01,dates, zRates_DV01);
    
    % Price the certificate with the shifted data 
    X_bucket = certificate_upfront(Certificate_data, dates, zRates_DV01,Data_capvol_DV01);
    
    % Compute the certificate bucket sensitivity 
    bucket_DV01(n_depos+ii) =  (X_bucket - X)/h * 1e-4; 

    % Compute the swap bucket sensitivity
    [bucket_swap_DV01(n_depos+ii,1),~, ~] = sensSwap(dates(1), swap_dates(1:2), swap_rates(1), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(n_depos+ii,2),~, ~] = sensSwap(dates(1), swap_dates(1:5), swap_rates(2), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(n_depos+ii,3),~, ~] = sensSwap(dates(1), swap_dates(1:10), swap_rates(3), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(n_depos+ii,4),~, ~] = sensSwap(dates(1), swap_dates(1:15), swap_rates(4), dates, discounts,discounts_DV01);
end

% Compute the bucket delta sensitivites for all swaps
for ii = 1:n_swaps

    % Increase the single mkt rate
    ratesSet_DV01 = ratesSet;
    ratesSet_DV01.swaps(ii,:) = ratesSet_DV01.swaps(ii,:) + h;


    % Do the bootstrap with new mkt rates
    [~, discounts_DV01] = bootstrap(datesSet, ratesSet_DV01);
    zRates_DV01 = zeroRates(dates,discounts_DV01)/100;

    % Compute the new Cap market prices and recalibrate the spot vol surface 
    Data_capvol_DV01 = Data_capvol;
    Data_capvol_DV01.cap_prices = Price_Cap_flat(Data_capvol_DV01, dates, zRates_DV01); 
    Data_capvol_DV01.sigma_spot = bootstap_vol(Data_capvol_DV01,dates, zRates_DV01);
    
    % Price the certificate with the shifted data 
    X_bucket = certificate_upfront(Certificate_data, dates, zRates_DV01,Data_capvol_DV01);
    
    % Compute the certificate bucket sensitivity 
    bucket_DV01(n_depos+n_futures+ ii) =  (X_bucket - X)/h * 1e-4;

    % Compute the swap bucket sensitivity
    [bucket_swap_DV01(n_depos+n_futures+ii,1),~, ~] = sensSwap(dates(1), swap_dates(1:2), swap_rates(1), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(n_depos+n_futures+ii,2),~, ~] = sensSwap(dates(1), swap_dates(1:5), swap_rates(2), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(n_depos+n_futures+ii,3),~, ~] = sensSwap(dates(1), swap_dates(1:10), swap_rates(3), dates, discounts,discounts_DV01);
    [bucket_swap_DV01(n_depos+n_futures+ii,4),~, ~] = sensSwap(dates(1), swap_dates(1:15), swap_rates(4), dates, discounts,discounts_DV01);
end
bucket_swap_DV01 = bucket_swap_DV01/h*1e-4;
t= toc;

fprintf(" Delta bucket sensitivities computed, elapsed time: %.2f \n\n",t);
fprintf("------------------------------------------------------------ \n\n");


%% Total vega


fprintf("------------------ Certificate total vega  ------------------\n\n");
% Set the increment for the derivative computation
h = 1e-8;

% Do the parallel shift of flat volatilities
Data_capvol_h = Data_capvol;
Data_capvol_h.flat_volatilities = Data_capvol_h.flat_volatilities +h;

% Compute the new Cap market prices and recalibrate the spot vol surface 
Data_capvol_h.cap_prices = Price_Cap_flat(Data_capvol_h, dates, zRates); 
Data_capvol_h.sigma_spot = bootstap_vol(Data_capvol_h,dates, zRates);

% Price the certificate with the shifted data
X_h_fwd = certificate_upfront(Certificate_data, dates, zRates,Data_capvol_h);

% Compute the certificate total vega 
total_vega = (X_h_fwd - X)/h *0.01;

% Compute also the total vega for 5y ATM Cap 
strike_1 = mean(ratesSet.swaps(5,:));
Expiry_date_1 = 5;

% Price the Cap both with shifted surface and mkt surface
Cap_price_1 = Cap_price_spot(Data_capvol, strike_1, Expiry_date_1, dates, discounts);
Cap_price_h_1 = Cap_price_spot(Data_capvol_h, strike_1, Expiry_date_1, dates, discounts);

% Compute the Cap total vega
total_vega_Cap_1 = (Cap_price_h_1-Cap_price_1)/h * 0.01;


fprintf(" The certificate total vega is: %.4f \n\n", total_vega);

%% Vega bucket sensitivities

% Import data of the Cap 
Expiry_date_2 = 15;
strike_2 = mean(ratesSet.swaps(13,:));
Cap_price_2 = Cap_price_spot(Data_capvol, strike_2, Expiry_date_2, dates, discounts);


% Initialize the vector of vega bucket for the certificate and 
% for the ATM Cap at 5y and 15y
vega_bucket = zeros(size(Data_capvol.flat_volatilities,1),1);
vega_bucket_cap = zeros(size(Data_capvol.flat_volatilities,1),2);

fprintf("------------------------------------------------------------ \n\n");

fprintf(" Starting computation of vega bucket sensitivities \n");

tic
% Compute the vega bucket
for ii = 1: size(Data_capvol.flat_volatilities,1)

    % Shift the flat volatilities of the ii_th expiry 
    Data_capvol_h = Data_capvol;
    Data_capvol_h.flat_volatilities(ii,:) = Data_capvol_h.flat_volatilities(ii,:) + h;
    
    % Compute the cap prices with the modified volatility surface
    Data_capvol_h.cap_prices = Price_Cap_flat(Data_capvol_h, dates, zRates); 
    
    % Recalibrate the spot volatility surface
    Data_capvol_h.sigma_spot = bootstap_vol(Data_capvol_h,dates, zRates);

    % Price the certificate with the new surface
    X_h_fwd = certificate_upfront(Certificate_data, dates, zRates,Data_capvol_h);
    
    % Compute the vega bucket for the certificate
    vega_bucket(ii) = (X_h_fwd - X)/h * 0.01;

    % Price the Caps with the new surface
    Cap_price_h_1 = Cap_price_spot(Data_capvol_h, strike_1, Expiry_date_1, dates, discounts);
    Cap_price_h_2 = Cap_price_spot(Data_capvol_h, strike_2, Expiry_date_2, dates, discounts);

    % Compute the vega bucket for the Caps
    vega_bucket_cap(ii, 1) = (Cap_price_h_1-Cap_price_1)/h * 0.01;
    vega_bucket_cap(ii, 2) = (Cap_price_h_2-Cap_price_2)/h * 0.01;

end
t= toc;
fprintf(" Vega bucket sensitivities computed, elapsed time: %.2f \n\n",t);
fprintf("------------------------------------------------------------ \n\n");


%% Coarse grained delta bucket


fprintf("-------- Hedge the certificate on the CG delta bucket --------\n\n");

% Initialize the vector of bucket dates
dates_bucket = [datesSet.depos(1:n_depos);datesSet.futures(1:n_futures,2); datesSet.swaps(1:n_swaps) ];

% Compute the weights of the first CGB
w_1 = interp1([datesSet.settlement,datesSet.swaps(2), datesSet.swaps(5)],[1,1,0],dates_bucket, "linear",0);

% Compute the weights of the second CGB 
w_2 = interp1([datesSet.swaps(2),datesSet.swaps(5), datesSet.swaps(10)],[0,1,0],dates_bucket,"linear",0);

% Compute the weights of the third CGB 
w_3 = interp1([datesSet.swaps(5),datesSet.swaps(10), datesSet.swaps(13)],[0,1,0],dates_bucket,"linear",0);

% Compute the weights of the fourth CGB
w_4 = interp1([datesSet.swaps(10),datesSet.swaps(13)],[0,1],dates_bucket,"linear",0);

% Initialize the matrix of weights
W = [w_1,w_2, w_3, w_4];

% Compute the Coarse grained delta bucket for the certificate
CGB_sens = W'*bucket_DV01 ;

% Compute the Coarse grained delta bucket for the swaps
CGB_sens_swap = W'*bucket_swap_DV01;

% Hedge the portfolio by solving the linear system
Hedging_portfolio_CGB_delta = -CGB_sens_swap\CGB_sens;

fprintf(" Position on 2y swap: %.2f € (%.2f)\n",Hedging_portfolio_CGB_delta(1)*Certificate_data.Notional,Hedging_portfolio_CGB_delta(1));
fprintf(" Position on 5y swap: %.2f € (%.2f)\n",Hedging_portfolio_CGB_delta(2)*Certificate_data.Notional,Hedging_portfolio_CGB_delta(2));
fprintf(" Position on 10y swap: %.2f € (%.2f)\n",Hedging_portfolio_CGB_delta(3)*Certificate_data.Notional,Hedging_portfolio_CGB_delta(3));
fprintf(" Position on 15y swap: %.2f € (%.2f)\n\n",Hedging_portfolio_CGB_delta(4)*Certificate_data.Notional,Hedging_portfolio_CGB_delta(4));



%% Hedge the certificate on total vega with a Cap 5y ATM

fprintf("----------- Hedge the certificate on the total vega -----------\n\n");

% Hedge the total vega of the certificate with the 5y ATM Cap
Hedging_portfolio_totalvega = - total_vega/total_vega_Cap_1;

fprintf(" Position on 5y ATM Cap: %.2f € (%.2f)\n\n",Hedging_portfolio_totalvega*Certificate_data.Notional, Hedging_portfolio_totalvega);

% Compute the total delta of the portfolio and hedge with a 5y swap

% Set the increment for the derivative computation
h = 1e-8;

% Shift all the mkt rates
ratesSet_DV01 = ratesSet;
ratesSet_DV01.depos = ratesSet_DV01.depos + h;
ratesSet_DV01.futures = ratesSet_DV01.futures + h;
ratesSet_DV01.swaps = ratesSet_DV01.swaps + h;

% Do the bootstrap
[~, discounts_DV01] = bootstrap(datesSet, ratesSet_DV01);
zRates_DV01 = zeroRates(dates,discounts_DV01)/100;

% Compute new cap prices and recalibrate spot volatility surface
Data_capvol_DV01.cap_prices = Price_Cap_flat(Data_capvol, dates, zRates_DV01); 
Data_capvol_DV01.sigma_spot = bootstap_vol(Data_capvol_DV01,dates, zRates_DV01);

% Price the certificate with new data
X_h = certificate_upfront(Certificate_data, dates, zRates_DV01,Data_capvol_DV01);

% Compute the total delta for portfolio assets
total_delta = (X_h - X)/h * 1e-4; 
[total_delta_swap_5y,~, ~] = sensSwap(dates(1), swap_dates(1:5), swap_rates(2), dates, discounts,discounts_DV01);
total_delta_swap_5y = total_delta_swap_5y/h *1e-4;
total_delta_Cap_1 = (Cap_price_spot(Data_capvol_DV01, strike_1, Expiry_date_1, dates, discounts_DV01) - Cap_price_1)/h *1e-4;

% Solve the linear system for the total Hedged portfolio 
Matrix_delta_vega = [total_vega_Cap_1, 0; total_delta_Cap_1, total_delta_swap_5y];
b_delta_vega = [total_vega; total_delta];

fprintf("-------- Hedge the certificate on the total delta-vega --------\n\n");

Hedging_portfolio_delta_vega = -Matrix_delta_vega\b_delta_vega;

fprintf(" Position on 5y ATM Cap: %.2f € (%.2f)\n",Hedging_portfolio_delta_vega(1)*Certificate_data.Notional,Hedging_portfolio_delta_vega(1));
fprintf(" Position on 5y swap: %.2f € (%.2f)\n\n",Hedging_portfolio_delta_vega(2)*Certificate_data.Notional,Hedging_portfolio_delta_vega(2));


%% Hedge the certificate on CGB vega with a Cap 5y ATM and a Cap 15y

fprintf("------------ Hedge the certificate on CGB vega ------------\n\n");
% Compute the weights of the first CGB
w_1 = interp1([Data_capvol.Settlement,Data_capvol.expiries(5),Data_capvol.expiries(12)],[1,1,0],Data_capvol.expiries, "linear",0);

% Compute the weights of the second CGB 
w_2 = interp1([Data_capvol.expiries(5),Data_capvol.expiries(12)],[0,1],Data_capvol.expiries,"linear",0);

% Initialize the matrix of weights
W = [w_1,w_2];

% Compute the CGB vega for Certificate and Caps
CGB_vega_sens_cap = W'*vega_bucket_cap;
CGB_vega_sens = W'*vega_bucket;

% Solve the linear system to Hedge the portfolio
Hedging_portfolio_CGB_vega = -CGB_vega_sens_cap\CGB_vega_sens;


fprintf(" Position on 5y ATM Cap: %.2f € (%.2f)\n",Hedging_portfolio_CGB_vega(1)*Certificate_data.Notional,Hedging_portfolio_CGB_vega(1));
fprintf(" Position on 15y ATM Cap: %.2f € (%.2f)\n\n",Hedging_portfolio_CGB_vega(2)*Certificate_data.Notional,Hedging_portfolio_CGB_vega(2));

%% Sensibility check
fprintf("------------------------------------------------------------\n\n")

fprintf(" Check on the sensibility computation:\n")
fprintf(" Total sens - sum of buckets\n\n");

fprintf(" Certificate sens \n")
fprintf(" delta: %.2d , %.2d\n", total_delta,sum(bucket_DV01) );
fprintf(" vega:  %.2d , %.2d \n\n", total_vega, sum(vega_bucket));

fprintf(" 5y swap sens\n")
fprintf(" delta: %.2d , %.2d\n\n", total_delta_swap_5y,sum(bucket_swap_DV01(:,2)) );

fprintf(" 5y ATM Cap sens\n")
fprintf(" vega: %.2d , %.2d\n", total_vega_Cap_1,sum(vega_bucket_cap(:,1)) );
















