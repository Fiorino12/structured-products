function price = Bachelier_LMM(B, delta_t,L_fwd, K,sigma, T)
% Compute the caplet price using the Bachelier model
%
% INPUT:
% B:        Discount factor
% delta_t:  time distance between fixing and payment   
% L_fwd:    Forward Libor
% K:        Strike rate
% sigma:    Libor volatility
% T:        Time to maturity
%
% OUTPUT:
% price:    caplet price

d_n = (L_fwd - K)./(sigma.*sqrt(T));

price = B.*delta_t.*((L_fwd - K).*normcdf(d_n) + sigma.*sqrt(T).*normpdf(d_n));
end