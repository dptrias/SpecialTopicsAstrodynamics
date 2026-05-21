function x_Lagrange = ClassicalLagrangePoints(mu, point)
    f_mu = @(x) x - (1 - mu) / (x + mu)^2 - mu / (x - (1 - mu))^2;
    df_mu = @(x) 1 - (1 - mu) * (1 / (x + mu)^3 - 3 * (x + mu) / (x + mu)^4) - mu * (1 / (x - (1 - mu))^3 - 3 * (x - (1 - mu)) / (x - (1 - mu))^4);

    if point == 1
        x0 = 0.5; % Initial guess for L1
    elseif point == 2
        x0 = 1.5; % Initial guess for L2
    else 
        fprintf('Invalid point. Please choose 1 for L1 or 2 for L2.\n');
        return;
    end

    % Newton-Raphson method
    tol = 1e-6;
    max_iter = 1000;
    x = x0;
    for i = 1:max_iter
        fx = f_mu(x);
        dfx = df_mu(x);
        if dfx == 0
            fprintf('Derivative is zero. No solution found.\n');
            return;
        end
        x_new = x - fx / dfx;
        if abs(x_new - x) < tol
            fprintf('L%d point found at x = %.6f\n', point, x_new);
            x_Lagrange = x_new;
            return;
        end
        x = x_new;
    end

end
