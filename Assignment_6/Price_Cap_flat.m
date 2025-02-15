function cap_prices = Price_Cap_flat(Data_capvol, dates, zRates)
% Price the Mkt Caps from flat volatilities
%
% INPUT: 
%  Data_capvol: struct with the following components:
%                   -Settlement:          Settlement date
%                   -strikes:             Row vector of strikes in the volatility matrix
%                   -expyear:             Column vector of expiries in the volatility matrix
%                                        expressed in yf (i.e. if the expiry is 1y then is reported 1)
%                   -expiries:            Column vector of the expiries in volatility
%                                        matrix in datenum
%                   -flat_volatilities:   mkt flat volatilities matrix
%  dates:       datenum dates of the bootstrap
%  zRates:      zRates obtained from the bootstrap
%
% OUTPUT:
%  cap_prices:  Matrix of the MKt Cap prices 

% Set the dates conventions
zRatesconvention = 3; % Act/365
capletconvention = 2; % Act/360

% Convert the bootstrap dates in year fraction with zRates convention
dates_yfz = yearfrac(dates(1), dates, zRatesconvention);

% Initialize the matrix of Cap prices
cap_prices = zeros(size(Data_capvol.flat_volatilities));

% Find the caplet payment dates for every Cap
caplet_dates= finddates(Data_capvol.Settlement, (3:3:Data_capvol.expyear(end)*12)',1);

% Convert the caplet dates with the convention used for pricing
caplet_dates_yfz = yearfrac(Data_capvol.Settlement, caplet_dates, zRatesconvention);
caplet_dates_yfc = yearfrac(Data_capvol.Settlement, caplet_dates, capletconvention);

% Compute the delta of the caplets
delta_caplet = diff(caplet_dates_yfc);

% Compute the discount factor in the caplet dates
disc_caplet = exp(-interp1(dates_yfz, zRates,caplet_dates_yfz).*caplet_dates_yfz);
    
% Compute the forward discount factors 
disc_fwd = disc_caplet(2:end)./disc_caplet(1:end-1);
    
% Compute the forward libor
L_fwd = 1./delta_caplet .*(1./disc_fwd-1);

% Compute the Mkt Cap prices
for ii=1:length(Data_capvol.expiries)
    
    % Find the upper index of caplets for the iith expiries
    u = Data_capvol.expyear(ii)*4;
    
    % Compute the cap prices summing the caplet prices
    cap_prices(ii,:) = sum(Bachelier_LMM(disc_caplet(2:u), delta_caplet(1:u-1),L_fwd(1:u-1), Data_capvol.strikes,Data_capvol.flat_volatilities(ii,:),caplet_dates_yfz(1:u-1)));

end
end


