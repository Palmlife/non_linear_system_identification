function [nb_opt, na_opt, AIC, whiteness, GIC] = model_order_selection(G_ML, f, fs, noise_var, na_max, nb_max, fmax, method)
% MODEL_ORDER_SELECTION - Select optimal model order using AICc, GIC, and whiteness test
%   method: 'LLS', 'TLS', 'GTLS', 'WLS', 'IWLS', 'IQML', or 'ML' (default: 'IQML')
%
%   Outputs:
%     nb_opt, na_opt: optimal orders (based on AICc)
%     AIC: AICc values for all (na, nb) combinations
%     whiteness: whiteness test results for all (na, nb) combinations  
%     GIC: Generalized Information Criterion values (more robust for colored noise)
%
%   AICc = N*log(V) + 2p + 2p(p+1)/(N-p-1)  -- assumes white Gaussian residuals
%   GIC  = N*log(V) + p*log(N)              -- BIC/MDL penalty, more conservative
%
%   GIC is preferred when residuals are not white (e.g., stochastic NL distortion)

    if nargin < 5, na_max = 6; end
    if nargin < 6, nb_max = 5; end
    if nargin < 7, fmax = fs/2; end
    if nargin < 8, method = 'IQML'; end
    
    idx = (f > 0) & (f < fmax);
    N = sum(idx); G = G_ML(idx); nv = noise_var(idx);
    
    AIC = nan(na_max, nb_max+1);
    GIC = nan(na_max, nb_max+1);
    whiteness = nan(na_max, nb_max+1);
    
    % Methods that don't need noise variance
    no_noise_methods = {'TLS', 'LLS', 'IWLS'};
    
    for na = 1:na_max
        for nb = 0:min(nb_max, na)
            try
                if any(strcmp(method, no_noise_methods))
                    [~, ~, Gm] = parametric_LS(G_ML, f, fs, nb, na, [], method, fmax);
                else
                    [~, ~, Gm] = parametric_LS(G_ML, f, fs, nb, na, noise_var, method, fmax);
                end
                res = (G - Gm(idx)) ./ sqrt(nv);
                V = mean(abs(res).^2);
                p = nb + 1 + na;  % number of parameters
                
                % AICc: Corrected Akaike Information Criterion
                AIC(na, nb+1) = aicc_criterion(V, p, N);
                
                % GIC: Generalized Information Criterion (BIC/MDL form)
                GIC(na, nb+1) = gic_criterion(V, p, N);
                
                whiteness(na, nb+1) = whiteness_test(res);
            catch
                % If method fails for this order, leave as NaN
                AIC(na, nb+1) = NaN;
                GIC(na, nb+1) = NaN;
                whiteness(na, nb+1) = NaN;
            end
        end
    end
    
    % Find optimal order (using AICc by default)
    [~, i] = min(AIC(:)); 
    [na_opt, nb_opt] = ind2sub(size(AIC), i); nb_opt = nb_opt - 1;
    
    % Also find GIC-optimal order for comparison
    [~, i_gic] = min(GIC(:));
    [na_gic, nb_gic] = ind2sub(size(GIC), i_gic); nb_gic = nb_gic - 1;
    
    fprintf('Optimal (%s): AICc -> na=%d, nb=%d (AICc=%.2f, Whiteness=%.0f%%)\n', ...
        method, na_opt, nb_opt, AIC(na_opt,nb_opt+1), whiteness(na_opt,nb_opt+1)*100);
    fprintf('              GIC  -> na=%d, nb=%d (GIC=%.2f, Whiteness=%.0f%%)\n', ...
        na_gic, nb_gic, GIC(na_gic,nb_gic+1), whiteness(na_gic,nb_gic+1)*100);
end

function w = whiteness_test(res)
% WHITENESS_TEST - Test if residuals behave like white noise
%   For white noise, autocorrelation at non-zero lags should be ~0.
%   Returns the fraction of ACF values within 95% confidence bounds.
%
%   w = 1.0 (100%): residuals appear white (good model)
%   w ~ 0.95: expected for truly white noise (5% false positives)
%   w < 0.5: residuals are correlated (model misspecified)

    % Stack real and imaginary parts into one real vector
    % (both parts should be white if the complex residual is white)
    r = [real(res); imag(res)]; 
    N = length(r);
    
    % Compute normalized autocorrelation function (ACF)
    % - Subtract mean to center the signal
    % - Compute up to min(20, N/4) lags (enough to detect correlation)
    % - 'normalized' scales so ACF(0) = 1
    [acf, lags] = xcorr(r-mean(r), min(20,floor(N/4)), 'normalized');
    
    % Test: for white noise with N samples, ACF at lag > 0 is approximately
    % N(0, 1/N), so 95% confidence interval is ±1.96/sqrt(N).
    % Count what fraction of positive-lag ACF values fall within bounds.
    w = mean(abs(acf(lags>0)) < 1.96/sqrt(N));
end

function g = gic_criterion(V, p, N)
% GIC_CRITERION - Generalized Information Criterion (BIC/MDL form)
%   Computes the GIC for model order selection.
%
%   Inputs:
%     V : Mean squared (normalized) residual, i.e., V = mean(|res|^2)
%     p : Number of model parameters (p = nb + 1 + na)
%     N : Number of data points (frequency bins)
%
%   Output:
%     g : GIC value (lower is better)
%
%   Formula:
%     GIC = N * log(V) + p * log(N)
%
%   Properties:
%     - More conservative than AICc (stronger penalty for complexity)
%     - Consistent: selects true order as N -> infinity
%     - More robust when residuals are colored (not white noise)
%     - Equivalent to BIC (Bayesian IC) and MDL (Minimum Description Length)
%
%   Comparison with AICc:
%     AICc penalty: 2p + 2p(p+1)/(N-p-1) ≈ 2p for large N
%     GIC penalty:  p * log(N)
%     For N > 8, log(N) > 2, so GIC penalizes complexity more strongly

    g = N * log(V) + p * log(N);
end

function a = aicc_criterion(V, p, N)
% AICC_CRITERION - Corrected Akaike Information Criterion
%   Computes the AICc for model order selection.
%
%   Inputs:
%     V : Mean squared (normalized) residual, i.e., V = mean(|res|^2)
%     p : Number of model parameters (p = nb + 1 + na)
%     N : Number of data points (frequency bins)
%
%   Output:
%     a : AICc value (lower is better)
%
%   Formula:
%     AICc = N * log(V) + 2p + 2p(p+1)/(N-p-1)
%
%   Properties:
%     - Correction term important when N/p is not large
%     - Assumes residuals are white Gaussian noise
%     - Tends to select higher orders than GIC/BIC
%     - Not consistent (may overfit as N -> infinity)

    a = N * log(V) + 2*p + 2*p*(p+1)/(N-p-1);
end