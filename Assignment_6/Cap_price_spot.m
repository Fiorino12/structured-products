function Cap_price = Cap_price_spot(Data_capvol, strike, Expiry_date, dates, discounts)
% Compute the Cap price with the spot volatilities with Bachelier formula
% 
% 
% INPUT:
%  Data_capvol:         struct with the following components:
%                           -Settlement:          Settlement date
%                           -strikes:             Row vector of strikes in the volatility matrix
%                           -expyear:             Column vector of expiries in the volatility matrix
%                                                 expressed in yf (i.e. if the expiry is 1y then is reported 1)
%                           -expiries:            Column vector of the expiries in volatility
%                                                 matrix in datenum
%                           -flat_volatilities:   mkt flat volatilities matrix
%                           -sigma_spot:          spot volatilities matrix
%  strike:              Cap strike price 
%  Expiry_date:         Cap expiry date 
%  dates:               datenum dates of the bootstrap
%  discounts:           discount factors obtained from the bootstrap
%
% OUTPUT:
%  Cap_price: price of the Cap


% Set the dates conventions
zRatesconvention = 3; % Act/365
capletconvention = 2; % Act/360

% Convert the bootstrap dates in year fraction with zRates convention
dates_yfz = yearfrac(dates(1), dates, zRatesconvention);

% Find the dates of expiry for Cap matrix
date_expiry_yfc = yearfrac(dates(1), finddates(dates(1), Data_capvol.expyear), capletconvention);

% Compute zRates
zRates = zeroRates(dates, discounts)/100;

% Find the caplet payment dates for the Cap
payment_dates = finddates(dates(1), (3:3:Expiry_date*12)',1);
payment_dates_yfz = yearfrac(dates(1), payment_dates, zRatesconvention);
payment_dates_yfc = yearfrac(dates(1), payment_dates, capletconvention);

% Compute the discounts for caplets
discounts_pay = exp(-interp1(dates_yfz, zRates, payment_dates_yfz).*payment_dates_yfz);

% Compute the caplet deltas
delta_t = yearfrac([dates(1); payment_dates(1:end-1)], payment_dates(1:end), capletconvention);

% Interpolate the sigma spot on caplet expiries and then on the strike
sigma_tot = interp1([0;date_expiry_yfc], [zeros(1,13);Data_capvol.sigma_spot], payment_dates_yfc(2:end) );
sigma_tot = interp1(Data_capvol.strikes', sigma_tot',strike ,'spline');

% Compute the forward Libor
disc_fwd = discounts_pay(2:end)./discounts_pay(1:end-1);
L_fwd = 1./delta_t(2:end) .*(1./disc_fwd-1);

% Compute the Cap price as the sum of caplets
Cap_price = sum(Bachelier_LMM(discounts_pay(2:end), delta_t(2:end),...
            L_fwd, strike,sigma_tot',  payment_dates_yfz(1:end-1)));



end