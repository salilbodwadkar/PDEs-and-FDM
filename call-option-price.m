% a MATLAB scheme to solve Black-Scholes using Crank-Nicholson

% S0 = Initial stock price
% K = strike price
% r = risk-free interest rate
% sigma = volatility
% t_mat = time to maturity
% T = number of time steps
% N = number of stock price steps

function optionPrice = crankNicolsonBS(S0, K, r, sigma, t_mat, T, N)
    dt = t_mat / T;
    dS = (2 * S0) / N;
    s = 0:dS:(2 * S0);
    t = 0:dt:t_mat;
  
    V = zeros(N+1, T+1);

    % Using standard european call BCs
    V(:, end) = max(s - K, 0);
    V(1, :) = 0;
    V(end, :) = 2 * S0 - K * exp(-r * t);

    % Crank-Nicolson coefficients
    alpha = 0.25 * dt * ((sigma^2 * s.^2) / dS^2 - r * s / dS);
    beta = -0.5 * dt * ((sigma^2 * s.^2) / dS^2 + r);
    gamma = 0.25 * dt * ((sigma^2 * s.^2) / dS^2 + r * s / dS);

    for j = T:-1:1
        % tri-diagonal
        A = diag(alpha(3:N-1), -1) + diag(1 - beta(2:N-1)) + diag(gamma(2:N-2), 1);
        
        b = V(2:N-1, j+1);
        b(1) = b(1) - alpha(2) * V(1, j);
        b(end) = b(end) - gamma(N-1) * V(N+1, j);

        V(2:N-1, j) = A \ b;
    end
    optionPrice = V;
end
