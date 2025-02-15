function [DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts,discounts_DV01)
% Compute the sensitivity indexes of an IRS (DV01, BPV, DVO1_z)
% 
% setDate:              settlement Date
% fixedLegPaymentDates: Payments Dates of the fixed Leg
% fixedRate:            swap rate
% dates:                dates of the discount factors computed in bootstrap
% discounts:            df obtained by the bootstrap
% discounts_DV01:       df of the curve by shifting the market date of 1bp

% Date for zero rates using act/365 convention

datemode = 3;
fixedLegdates_yf_3= yearfrac(setDate, fixedLegPaymentDates, datemode);
dates_yf = yearfrac(setDate, dates, datemode);

% Compute the zRates corresponding to the fixed leg 
% for both shifted curve and real curve

zRates= zeroRates(dates, discounts);
zRates_DV01= zeroRates(dates, discounts_DV01);

% Compute zero rates, then DF for both curves

zRates_fixedLeg=interp1(dates_yf,zRates./100,fixedLegdates_yf_3);
zRates_fixedLeg_DV01=interp1(dates_yf,zRates_DV01./100,fixedLegdates_yf_3);

discounts_fl= exp(-fixedLegdates_yf_3.*zRates_fixedLeg);
discounts_fl_DV01= exp(-fixedLegdates_yf_3.*zRates_fixedLeg_DV01);

% Compute delta_t then BPV for both curves and 
% date for zero rates using act/365 convention

datemode = 6;

delta_t= yearfrac([setDate;fixedLegPaymentDates(1:end-1)], fixedLegPaymentDates, datemode);
BPV=sum(delta_t.*discounts_fl);
BPV_DV01=sum(delta_t.*discounts_fl_DV01);

% Compute NPV for both curves and find the difference (DV01)

NPV = BPV*fixedRate/100-1+discounts_fl(end);
NPV_DV01 = BPV_DV01*(fixedRate)/100 - 1 + discounts_fl_DV01(end);
DV01 = NPV_DV01-NPV; 

% compute DV01_z by shifting the zero rates

zRates_fixedLeg_DV01z=zRates_fixedLeg+0.0001;
discounts_fl_DV01z= exp(-fixedLegdates_yf_3.*zRates_fixedLeg_DV01z);
BPV_DV01z=sum(delta_t.*discounts_fl_DV01z);
NPV_DV01z= BPV_DV01z*fixedRate/100 - 1 + discounts_fl_DV01z(end);
DV01_z = NPV_DV01z-NPV;










end