% Define parameter

n = 10;
Hm = rand(1) + .5;
tau = rand(n, 1) + 3;
B = zeros(n);
A = zeros(n, n, n);
f = rand(n) + 1.5 + diag(10 * rand(n, 1) + .5);
for k = 1:n
    B(:, k) = rand(n, 1);
    B(k, k) = 0;
    B(:, k) = B(:, k) / sum(B(:, k));
    A(:, :, k) = diag(f(:, k) - Hm);
    A(k, :, k) = ones(1, n) * Hm;
    A(k, k, k) = f(k, k);
end
[V, D] = eig(B);
p = V(:, 1) / sum(V(:, 1));
T = 1000;
x0 = ones(n, 1) / n;
nT = 10;


% Do simulation

[tx, xn, logN, xi] = SimulatePopulationDynamics(x0, T, A, B, tau, nT);


% Plot result and compare to theory

close all;
plot(tx, logN, 'b');
hold on;
DeltaR = repmat(diag(f)', n, 1) - f;
gammaR = (sum(p .* tau .* diag(f)) - sum(B .* log(1 + DeltaR' / Hm) * p)) / sum( p .* tau);
plot(tx, tx * gammaR, 'r');


% Save data

% csvwrite('SimRespSw.csv', [tx', logN']);
