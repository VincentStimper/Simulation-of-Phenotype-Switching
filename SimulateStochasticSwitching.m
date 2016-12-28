% Define parameter

n = 10;
H = (.5 * rand(n) + .5) / n / 5;
H = H - diag(sum(H));
tau = rand(n, 1) + 3;
B = zeros(n);
A = zeros(n, n, n);
f = rand(n) + 0.5 + diag(10 * rand(n, 1) + .5);
for k = 1:n
    B(:, k) = rand(n, 1);
    B(k, k) = 0;
    B(:, k) = B(:, k) / sum(B(:, k));
    A(:, :, k) = diag(f(:, k)) + H;
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
Delta = repmat(diag(f)', n, 1) - f;
DeltaS = 1 ./ (1 ./ Delta + 1 ./ Delta');
gammaS = (sum(p .* tau .* diag(f)) - sum(p .* tau .* diag(H)) - sum(B .* log(1 + DeltaS ./ H) * p)) / sum( p .* tau);
plot(tx, tx * gammaS, 'r');


% Save data

% csvwrite('SimStochSw.csv', [tx', logN']);
