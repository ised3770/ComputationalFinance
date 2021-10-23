% Option Pricing with Explicit Scheme
clear
clc
% -------- Parameters -------- %
K = 15;
r = 0.1;
sigma = 0.25;
T = 0.5; % Final time
gamma = 1;
S_max = 4*K;
N = 1000;  % 140
M = 100;  % 50
dt = T/N;     % Time step
ds = S_max/M; % Price step
j = 0:M;      % Stock (price) points
%j = 1:M-1
n = 0:N;      % Time points

% Initialize grid
V(1:M+1, 1:N+1) = 0;

% Coeffieicents of tridiagonal matrix A
j2 = j.*j;
sigma2 = sigma*sigma;
j_gam = (j.*ds).^(2*gamma);
% Note S_j = j.*ds (price vector)

l_j = 0.5*sigma2*j_gam*(dt/(ds)^2) - j.*r*0.5*dt;
d_j = 1 - sigma2 * j_gam*(dt/(ds)^2) - dt*r;
u_j = 0.5*sigma2*j_gam*(dt/(ds)^2) + j.*r*0.5*dt;

% Initial condition, payoff function at expiry t = T
% V(j, N) = max(S_j-K ,0)
V(:, end) = max((j.*ds)-K, 0);

time = (0:N)*dt;
% Boundrary conditions
V(1, :) = 0;
%V(end, :) = S_max-K*exp(-r*flip(n*dt)); % flip?
V(end, :) = S_max-K*exp(-r*(T-time));


% Tridiagonal matrix (size: M-1 x N-1)
A =  diag(l_j(3:M), -1) + diag(d_j(2:M)) + diag( u_j(2:M-1), 1);

% ----- Main Loop ------ %
for n = flip(1:1:N) % Decreasing order
    % V_n-1 = A*V_n + b
    V(2:end-1, n) = A*V(2:end-1, n+1);
    % Add remaining terms lj(1), uj(end) from b vector
    V([2 end-1], n) = V([2 end-1], n) + [l_j(2) u_j(end)]' .* V([1 end], n+1);
    %V([2 end-1], n) = V([2 end-1], n) + [l_j(1) u_j(end)]' .* V([1 end], n+1);
end

time = (0:N)*dt;
S = j.*ds; % Price vector

% ----- Exact solution ------ %
V_exact = bsexact_vec(sigma, r, K, time, S);

% ----- Plots ------ %
figure(1)
plot(S, V(:, 1), 'LineWidth',2)
hold on
plot(S, max((j.*ds)-K, 0))
hold off
legend('Option Value/price','Payoff function') 
% 
figure(2)
plot(S, V(:, 1), 'r-') % Red: Value when contract is signed t=0
hold on
plot(S, V(:, round(N/2)), 'g-') % Green: Value at half maturit T/2
hold on
plot(S, V(:, N+1), 'b-' ) % Blue: Value at expiry, matury t=T
hold on
plot(S, V_exact(:, N+1), '*-') % Exact sol.
xlabel('S')
ylabel('V(s,t)')
xlim([0 20])
legend('Value at t=0', 'Value at t=T/2', 'Value at t=T (expiry)', 'Exact sol.')
% 
% figure(3)
% mesh(time, S, V)