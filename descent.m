function [opt_x, opt_f, opt_df, opt_alpha] = descent(dimension)
format compact

tol = 1e-6;
maxIter = 1000;
maxLs = 100;

opt_f = inf(1, 1);

for count = 1:50,
    k = 1;
    x(:, 1) = randn(dimension, 1);
    [f(k), df(:, k)] = obj(x(:, k));

    while norm(df(:, k)) > tol && k <= maxIter,
        [x(:, k + 1), alpha(k)] = inexact(x(:, k), f(k), df(:, k), maxLs);
        k = k + 1;
        [f(k), df(:, k)] = obj(x(:, k));
    end
    if opt_f(end) > f(end),
        opt_x = x;
        opt_f = f;
        opt_df = df;
        opt_alpha = alpha;
    end
end

if dimension == 2,
    draw_levy1;
    plot(opt_x(1, :), opt_x(2, :), 'r-');
end


function [x, alphaK] = inexact(xK, fK, dfK, maxLs)
    ds = -dfK;
    alphaK = 1;
    x = xK + alphaK * ds;
    f = obj(x);
    l = 1;
    while f >= fK && l <= maxLs,
        alphaK = alphaK / 2;
        x = xK + alphaK * ds;
        f = obj(x);
        l = l + 1;
    end
    if l > maxLs,
        error('Line Search Error');
    end
    [f, df] = obj(x);
    