function [tx, xn, logN, xi] = SimulatePopulationDynamics(x0, T, A, B, tau, nT)
    % Simulation of Population Dynamics acording to the model taken from
    % E. Kussell and S. Leibler. Phenotypic diversity, population growth,
    % and information in fluctuating environments. Science,
    % 309(5743):2075–2078, 2005.
    %
    % Input
    % A     A(:, :, k) is matrix A_k in paper
    % B     B(i, j) are probabilities b_{ij} in paper
    % tau   tau(i) is t_i in paper, interval lengths exponentially
    %       distributed
    % T     time length of simulation
    % x0    column vector with initial values of population vector
    % 
    % Output
    % tx    time
    % xn    xn(:,i) is population vector at time t(i) divided by population
    %       size
    % logN  logN(i) is logarithm of population size at time t(i)
    % xi    xi(i) index of environment at time t(i)
    % nT    number of time points simulated in one environment
    
    tx = 0;
    xn = x0;
    n = length(tau);
    xi(1) = ceil(rand(1) * n);
    N = sum(x0);
    logN = log(N);
    for k = 1:n
        pdEnv(k) = makedist('Multinomial', 'Probabilities', B(:, k));
    end
    
    while tx(end) < T
        dt = exprnd(tau(xi(end)));
        if tx(end) + dt < T
            tx = [tx, tx(end) + dt / nT * (1 : nT)];
        else
            tx = [tx, tx(end) + (T - tx(end)) / nT * (1 : nT)];
        end
        M = expm(A(:, :, xi(end)) * (tx(end) - tx(end - 1)));
        for k = 1 : nT
            xTil = M * xn(:, end);
            s = sum(xTil);
            logN = [logN, logN(end) + log(s)];
            xn = [xn, xTil / s]; 
        end
        xi = [xi, random(pdEnv(xi(end)), 1)];
    end

end

