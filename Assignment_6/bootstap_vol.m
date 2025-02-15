function sigma_spot = bootstap_vol(Data_capvol,dates, zRates)
% Bootstrap the spot volatilities from mkt data
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
%                           -cap_prices:          matrix of mkt cap prices
%  dates:       datenum dates of the bootstrap
%  zRates:      zRates obtained from the bootstrap
%
% OUTPUT:
%  sigma_spot: spot volatilities matrix 

% Set the dates conventions
zRatesconvention = 3; % Act/365
capletconvention = 2; % Act/360

% Convert the bootstrap dates in year fraction with zRates convention
dates_yfz = yearfrac(dates(1), dates, zRatesconvention);

% Compute the yf of spot volatilities expiry dates
cap_dates_yfc = yearfrac(Data_capvol.Settlement, Data_capvol.expiries, capletconvention);

% Compute the difference of Mkt Cap prices
delta_cap_price = diff(Data_capvol.cap_prices);

% Initialize the matrix of spot volatilities
n = size(delta_cap_price,1)+1;
m = size(delta_cap_price,2);
sigma_spot = zeros(n, m);

% Spot volatilities of the first year are the same of the flat ones
sigma_spot(1,:) = Data_capvol.flat_volatilities(1,:);

% Find the dates of all caplets and convert in yf with different
% conventions
caplet_dates= finddates(Data_capvol.Settlement, (Data_capvol.expyear(1)*12:3:Data_capvol.expyear(end)*12)',1);
caplet_dates_yfz = yearfrac(Data_capvol.Settlement, caplet_dates, zRatesconvention);
caplet_dates_yfc = yearfrac(Data_capvol.Settlement, caplet_dates, capletconvention);

% Compute the delta of the caplets
delta_caplet = diff(caplet_dates_yfc);

% Compute the caplets discount factors 
disc_caplet = exp(-interp1(dates_yfz, zRates,caplet_dates_yfz).*caplet_dates_yfz);

% Compute the forward Libor
disc_fwd = disc_caplet(2:end)./disc_caplet(1:end-1);
L_fwd = 1./delta_caplet .*(1./disc_fwd-1);

% Compute the spot vol 
for ii = 2:n 
   
   % Find the indexes of the caplets in the delta cap
   l= Data_capvol.expyear(ii-1)*4-2;
   u = Data_capvol.expyear(ii)*4-3;
    
   % Solve the system for mkt delta cap == cap with spot vol
   sigma_f = @(s) [sigma_spot(1:ii-1,:); s] ;
   sigma_cap = @(s) interp1(cap_dates_yfc(1:ii), sigma_f(s), caplet_dates_yfc(l:u));
   delta_cap_price_spot = @(s) sum(Bachelier_LMM(disc_caplet(l:u), delta_caplet(l-1:u-1),...
   L_fwd(l-1:u-1), Data_capvol.strikes,sigma_cap(s),  caplet_dates_yfz(l-1:u-1)));
   sigma_spot(ii,:) = fsolve(@(s) delta_cap_price_spot(s) - delta_cap_price(ii-1,:),sigma_spot(ii-1,:), optimoptions('fsolve','Display','off','OptimalityTolerance',1e-8) ); 
  
end

end