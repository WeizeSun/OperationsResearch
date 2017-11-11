function [opt_x, opt_f, opt_k, opt_l, opt_beta] = newton(dimension)
format compact

tol = 1e-6;
maxIter = 100;
maxLs = 10;

opt_f = inf(1,1);

for count = 1:100,
    x(:, 1) = randn(dimension, 1);

    k = 1;
    [f(k), df(:, k), ddf(:, :, k)] = obj(x(:, k));

    while norm(df(:, k)) >= tol && k <= maxIter,

        beta(k) = 0;
        l = min(eig(ddf(:,:,k)));
        if l < 0,
            beta(k) = floor(1 - l);
        end
    
        dN = -inv(ddf(:, :, k) + beta(k) * eye(length(x(:,k)))) * df(:, k);
        x(:, k+1) = dN + x(:, k);
        f(k+1) = obj(x(:,k+1));
    
        l(k) = 1;
        while l(k) < maxLs && f(k+1) >= f(k),
            l(k) = l(k) + 1;
            x(:, k+1) = x(:, k) + 2^-l(k) * dN;
            f(k+1) = obj(x(:,k+1));
        end
    
        if l(k) > maxLs,
            error('Line search error');
        end
    
        k = k + 1;
        [f(k), df(:, k), ddf(:, :, k)] = obj(x(:,k));

    end
    if opt_f(end) > f(end),
        opt_f = f;
        opt_x = x;
        opt_l = l;
        opt_beta = beta;
        opt_k = k;
    end
end

if dimension == 2,
    draw_levy1;
    plot(opt_x(1,:) ,opt_x(2,:), 'r-o');
end
