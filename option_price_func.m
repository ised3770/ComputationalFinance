clc
clear all

% ------ Variables ------ %
S0 = 14;
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5;
gamma = 1;
delta_t = 0.0001; % Time step size
N = 100; % Number of simulations

% ------ Initialize price and time ------ %
V = zeros(1, N); % Vector with value/price of option
S_vec = zeros(1, N);

% ------ Loop ------ %
for i = 1:N
    t = 0; 
    S = S0;
    while t < T
        t = t + delta_t;
        S = S + r * S * delta_t + sigma * (S^gamma) * randn * sqrt(delta_t);
    end
    S_vec(i) = S;
    V(i) = max(S - K, 0);
end

price_at_0 = exp(-r*T)*mean(V); % Discount
disp(price_at_0)
%plot([1:N], S_vec)

