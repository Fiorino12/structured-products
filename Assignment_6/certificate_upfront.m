function X = certificate_upfront(Certificate_data, dates, zRates,Data_capvol)
% Compute the upfront(%) of a certificate on the Libor rates 
% 
% INPUT: 
%  Certificate_data:    struct with the following components:
%                           -Notional:     Notional of the certificate
%                           -T:            Expiry of the certificate (in years)
%                           -spol_A:       spread over libor of part A
%                           -spol_B:       spread over libor of part B
%                           -first_libor:  first Libor fixed when the
%                                          contract is done
%                           -strike_5:     threshold of the part B capped
%                                          libor (from 0 to 5y)
%                           -strike_5_10:  threshold of the part B capped
%                                          libor (from 5 to 10y)
%                           -strike_10_15: threshold of the part B capped
%                                          libor (from 0 to 5y)
%  dates:               datenum dates of the bootstrap
%  zRates:              zRates obtained from the bootstrap
%  Data_capvol:         struct with the following components:
%                           -Settlement:          Settlement date
%                           -strikes:             Row vector of strikes in the volatility matrix
%                           -expyear:             Column vector of expiries in the volatility matrix
%                                                 expressed in yf (i.e. if the expiry is 1y then is reported 1)
%                           -expiries:            Column vector of the expiries in volatility
%                                                 matrix in datenum
%                           -flat_volatilities:   mkt flat volatilities matrix
%                           -sigma_spot:          spot volatilities matrix
%
% OUTPUT:
%  cap_prices:  Matrix of the MKt Cap prices 

% Set the dates conventions
zRatesconvention = 3; % Act/365
capletconvention = 2; % Act/360

% Initialize the array of Part B strikes 
strikes = [ones(1,19)*Certificate_data.strike_5, ones(1,20)*Certificate_data.strike_5_10, ones(1,20)*Certificate_data.strike_10_15];

% Convert the bootstrap dates in year fraction with zRates convention
dates_yfz = yearfrac(dates(1), dates, zRatesconvention);

% Find the dates of expiry for Cap matrix
date_expiry_yfc = yearfrac(dates(1), finddates(dates(1), Data_capvol.expyear), capletconvention);

% Find caplet payment dates
payment_dates = finddates(dates(1), (3:3:Certificate_data.T*12)',1);
payment_dates_yfz = yearfrac(dates(1), payment_dates, zRatesconvention);
payment_dates_yfc = yearfrac(dates(1), payment_dates, capletconvention);

% Find caplet discounts
discounts_pay = exp(-interp1(dates_yfz, zRates, payment_dates_yfz).*payment_dates_yfz);

% Compute the libor deltas
delta_t = yearfrac([dates(1); payment_dates(1:end-1)], payment_dates(1:end), capletconvention);

% Compute the BPV for part A
BPV = sum(delta_t.*discounts_pay);

% Compute the part A NPV
NPV_A = 1-discounts_pay(end) + Certificate_data.spol_A*BPV;

% Interpolate the sigma spot on caplet expiries and then on strikes
sigma_tot = interp1([0;date_expiry_yfc], [zeros(1,13);Data_capvol.sigma_spot], payment_dates_yfc );
sigma_tot = interp1(Data_capvol.strikes', sigma_tot', strikes'-Certificate_data.spol_B, 'spline');

% Select the volatility
sigma_tot = [sigma_tot(1, 2:20),sigma_tot(2, 21:40), sigma_tot(3, 41:60) ];

% Compute the fwd discounts and fwd libors
disc_fwd = discounts_pay(2:end)./discounts_pay(1:end-1);
L_fwd = 1./delta_t(2:end) .*(1./disc_fwd-1);

% Compute the caplet NPVs of Part B
Caplet_NPV= sum(Bachelier_LMM(discounts_pay(2:end), delta_t(2:end),...
            L_fwd, strikes'-Certificate_data.spol_B,sigma_tot',  payment_dates_yfz(1:end-1)));

% Compute Part B NPVs
NPV_B = Certificate_data.first_libor * delta_t(1)*discounts_pay(1) + ...
    Certificate_data.spol_B * sum(delta_t(2:end).*discounts_pay(2:end)) + ...
    discounts_pay(1) - discounts_pay(end) - Caplet_NPV;

% Compute the upfront value
X = NPV_A-NPV_B;


end