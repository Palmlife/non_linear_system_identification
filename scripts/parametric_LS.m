function [B, A, G_model, n_iter] = parametric_LS(G_ML, f, fs, nb, na, noise_var, method, fmax)
% Parametric TF estimation: G(z) = B(z)/A(z)
% Methods: 'LLS', 'WLS', 'IWLS', 'IQML' (monic constraint: a0=1)
%          'TLS' (norm constraint: ||theta||=1)
%          'GTLS' (Generalized TLS: ||C_J * theta|| = 1, uses GSVD)
%          'NLS' or 'ML' (Nonlinear LS / Maximum Likelihood via Levenberg-Marquardt)
%
% IWLS: Iterative WLS with weights 1/|A(z)|^2 (Sanathanan-Koerner, no noise info)
% IQML: Iterative QML with weights 1/(sigma^2 * |A(z)|^2) (consistent, needs noise)
%
% n_iter: number of iterations (0 for non-iterative methods)

    if nargin < 7, method = []; end
    if nargin < 6, noise_var = []; end
    if nargin < 8 || isempty(fmax), fmax = fs/2; end
    if isempty(method)
        if ~isempty(noise_var), method = 'WLS'; 
        else, method = 'LLS'; end
    end
    
    % Select positive frequencies up to fmax
    idx = (f > 0) & (f < fmax);
    G = G_ML(idx);
    z = exp(1j*2*pi*f(idx)/fs);
    w = ones(size(G)); 
    if ~isempty(noise_var), w = 1./(abs(noise_var(idx)) + eps); end
    
    n_iter = 0;  % Default for non-iterative methods
    
    switch upper(method)
        case 'TLS'
            % Total LS: constraint ||theta||=1, use SVD
            Phi_tls = [z.^(-(0:nb)), -G.*z.^(-(0:na))];
            Phi_r = [real(Phi_tls); imag(Phi_tls)];
            [~,~,V] = svd(Phi_r, 'econ');
            theta_full = V(:,end);  % Right singular vector for smallest singular value
            B = theta_full(1:nb+1).';
            A = theta_full(nb+2:end).';
            % Normalize so A(1) = 1 for consistent output
            B = B / A(1); 
            A = A / A(1);
            
        case 'GTLS'
            % Generalized Total Least Squares using GSVD
            % Minimizes ||J * theta||^2 subject to ||C_J^{1/2} * theta||^2 = 1
            % Uses Generalized SVD of (J, C_J^{1/2})
            
            Z_b = z.^(-(0:nb));           % Numerator regressors
            Z_a = z.^(-(0:na));           % Denominator frequency terms
            
            % Build regressor matrix J (homogeneous form)
            J = [Z_b, -G.*Z_a];
            
            % Number of parameters
            n_b = nb + 1;
            n_a = na + 1;
            n_theta = n_b + n_a;
            K = length(G);
            
            % Get sqrt of noise standard deviation for C_J^{1/2}
            % noise_var = sigma^2, so C_J^{1/2} uses sigma^{1/2} = (noise_var)^{1/4}
            if ~isempty(noise_var)
                sqrt_sigma = abs(noise_var(idx)).^(1/4);
            else
                sqrt_sigma = ones(K, 1);
            end
            
            % Build C_J^{1/2}: square root of noise covariance structure matrix
            % The noise only affects columns containing G (denominator part)
            % C_J^{1/2}(k,:) = sigma_k^{1/2} * [0, ..., 0, z_k^0, z_k^{-1}, ..., z_k^{-na}]
            sqrt_C_J = zeros(K, n_theta);
            sqrt_C_J(:, n_b+1:end) = diag(sqrt_sigma) * Z_a;
            
            % Convert to real representation for GSVD
            J_r = [real(J); imag(J)];
            sqrt_C_J_r = [real(sqrt_C_J); imag(sqrt_C_J)];
            
            % Generalized SVD: J_r = U * C * X' and sqrt_C_J_r = V * S * X'
            % MATLAB gsvd: [U,V,X,C,S] = gsvd(A,B) gives A = U*C*X' and B = V*S*X'
            % Generalized singular values: gamma_i = C(i,i) / S(i,i)
            % Solution: column of X^{-T} corresponding to smallest gamma
            [U, V, X, C_gsvd, S_gsvd] = gsvd(J_r, sqrt_C_J_r);
            
            % Extract diagonal elements
            % C_gsvd and S_gsvd may have different sizes
            [m_c, n_c] = size(C_gsvd);
            [m_s, n_s] = size(S_gsvd);
            
            % Compute generalized singular values
            gamma = inf(n_theta, 1);
            for i = 1:min([m_c, n_c, m_s, n_s, n_theta])
                if abs(S_gsvd(i,i)) > 1e-14
                    gamma(i) = abs(C_gsvd(i,i)) / abs(S_gsvd(i,i));
                elseif abs(C_gsvd(i,i)) < 1e-14
                    gamma(i) = 0;  % Both zero - degenerate case
                end
            end
            
            % Find the smallest generalized singular value
            [~, min_idx] = min(gamma);
            
            % Solution is X^{-T}(:,i) / s_i to satisfy ||C_J^{1/2} * theta||^2 = 1
            % Since J = U*C*X' and C_J^{1/2} = V*S*X', we need theta = inv(X')(:,i) / s_i
            X_inv_T = inv(X');
            s_i = S_gsvd(min_idx, min_idx);
            theta_full = real(X_inv_T(:, min_idx)) / s_i;
            
            B = theta_full(1:n_b).';
            A = theta_full(n_b+1:end).';
            
            % Normalize so A(1) = 1 for consistent output
            B = B / A(1);
            A = A / A(1);
            
        case 'LLS'
            Phi = [z.^(-(0:nb)), -G.*z.^(-(1:na))];
            theta = solve_ls(Phi, G, ones(size(G)));
            B = theta(1:nb+1).';
            A = [1, theta(nb+2:end).'];
            
        case 'WLS'
            Phi = [z.^(-(0:nb)), -G.*z.^(-(1:na))];
            theta = solve_ls(Phi, G, w);
            B = theta(1:nb+1).';
            A = [1, theta(nb+2:end).'];
            
        case 'IWLS'
            % Iterative Weighted Least Squares (Sanathanan-Koerner)
            % Weights by 1/|A(z)|^2 only (no noise variance)
            % Not consistent, but doesn't require noise information
            Phi = [z.^(-(0:nb)), -G.*z.^(-(1:na))];
            theta = solve_ls(Phi, G, ones(size(G)));
            for iter = 1:100
                A = [1, theta(nb+2:end).'];
                A_z = z.^(-(0:na)) * A.';
                w_iwls = 1 ./ (abs(A_z).^2 + eps);  % No noise weighting
                theta_new = solve_ls(Phi, G, w_iwls);
                rel_change = norm(theta_new - theta) / (norm(theta) + eps);
                fprintf('IWLS iter %d: rel_change = %.2e\n', iter, rel_change);
                if rel_change < 1e-8, break; end
                theta = theta_new;
            end
            n_iter = iter;
            B = theta(1:nb+1).';
            A = [1, theta(nb+2:end).'];
            
        case 'IQML'
            % Iterative Quadratic ML - weights by 1/(sigma^{2r} * |A(z)|^{2r})
            % Consistent estimator (requires noise variance)
            % Parameter r in [0,1]: r=1 is ML, r=0 is LLS
            r_iqml = 1;  % Set r=1 for ML estimate (can be adjusted)
            Phi = [z.^(-(0:nb)), -G.*z.^(-(1:na))];
            theta = solve_ls(Phi, G, ones(size(G)));
            for iter = 1:100
                A = [1, theta(nb+2:end).'];
                A_z = z.^(-(0:na)) * A.';
                % Equation error variance: sigma_e^2 = sigma^2 * |A(z)|^2
                % Weight: 1 / (sigma_e^2)^r = 1 / (sigma^{2r} * |A(z)|^{2r})
                sigma_e_sq = (1./w) .* (abs(A_z).^2);  % w = 1/sigma^2
                w_iqml = 1 ./ (sigma_e_sq.^r_iqml + eps);
                theta_new = solve_ls(Phi, G, w_iqml);
                rel_change = norm(theta_new - theta) / (norm(theta) + eps);
                fprintf('IQML iter %d: rel_change = %.2e\n', iter, rel_change);
                if rel_change < 1e-8, break; end
                theta = theta_new;
            end
            n_iter = iter;
            B = theta(1:nb+1).';
            A = [1, theta(nb+2:end).'];
                    
        case {'NLS', 'ML'}
            % Nonlinear Least Squares via Levenberg-Marquardt
            % Minimizes weighted output error: ||W*(G - B(z)/A(z))||^2
            % where W = diag(1/sigma_k) weights by inverse noise std
            % Uses full parameterization with ||theta||=1 constraint
            % theta = [b0, b1, ..., bnb, a0, a1, ..., ana]
            
            % Initialize with LLS (convert to full parameterization)
            Phi = [z.^(-(0:nb)), -G.*z.^(-(1:na))];
            theta_lls = solve_ls(Phi, G, ones(size(G)));
            % theta_lls = [b0,...,bnb, a1,...,ana] (monic: a0=1 implicit)
            % Convert to full: [b0,...,bnb, a0=1, a1,...,ana]
            theta = [theta_lls(1:nb+1); 1; theta_lls(nb+2:end)];
            theta = theta / norm(theta);  % Normalize to ||theta||=1
            
            % Levenberg-Marquardt iterations
            K = length(G);
            Z_b = z.^(-(0:nb));    % [K x (nb+1)]
            Z_a = z.^(-(0:na));    % [K x (na+1)] - full, including a0
            n_params = (nb + 1) + (na + 1);  % All b's and all a's
                    
                    % Weight vector: 1/sigma (inverse of noise std)
                    if ~isempty(noise_var)
                        W_vec = 1 ./ sqrt(abs(noise_var(idx)) + eps);
                    else
                        W_vec = ones(K, 1);
                    end
                    
                    % LM parameters
                    lambda_up = 10;
                    lambda_down = 5;
                    lambda_min = 1e-12;
                    lambda_max = 1e12;
                    
                    % Compute initial cost (weighted)
                    B_curr = theta(1:nb+1);
                    A_curr = theta(nb+2:end);
                    B_z = Z_b * B_curr;
                    A_z = Z_a * A_curr;
                    G_model_curr = B_z ./ A_z;
                    epsilon = G - G_model_curr;
                    epsilon_w = W_vec .* epsilon;  % Weighted error
                    cost_old = sum(abs(epsilon_w).^2);
                    cost_init = cost_old;
                    
                    % Compute initial Jacobian for lambda initialization
                    J_b = -Z_b ./ A_z;
                    J_a = G_model_curr .* Z_a ./ A_z;
                    J = [J_b, J_a];
                    J_w = W_vec .* J;  % Weighted Jacobian
                    J_real = [real(J_w); imag(J_w)];
                    lambda = 1e-3 * svds(J_real, 1);
                    
                    n_lm_iter = 0;
                    best_theta = theta;
                    best_cost = cost_old;
                    
                    for iter = 1:100
                        % Extract current B and A coefficients
                        B_curr = theta(1:nb+1);
                        A_curr = theta(nb+2:end);
                        
                        % Evaluate current model: G_model = B(z)/A(z)
                        B_z = Z_b * B_curr;
                        A_z = Z_a * A_curr;
                        G_model_curr = B_z ./ A_z;
                        
                        % Weighted output error: epsilon_w = W * (G - G_model)
                        epsilon = G - G_model_curr;
                        epsilon_w = W_vec .* epsilon;
                        
                        % Weighted Jacobian: J_w = W * J
                        J_b = -Z_b ./ A_z;                   % [K x (nb+1)]
                        J_a = G_model_curr .* Z_a ./ A_z;    % [K x (na+1)]
                        J = [J_b, J_a];                      % [K x n_params]
                        J_w = W_vec .* J;                    % Weighted Jacobian
                        
                        % Convert to real (stack real and imaginary parts)
                        epsilon_real = [real(epsilon_w); imag(epsilon_w)];
                        J_real = [real(J_w); imag(J_w)];
                        
                        % SVD of weighted Jacobian: J_w = U * Sigma * V'
                        [U_svd, Sigma_svd, V_svd] = svd(J_real, 'econ');
                        sigma_vec = diag(Sigma_svd);
                        
                        % Store cost before this iteration for convergence check
                        cost_before_iter = cost_old;
                        
                        step_accepted = false;
                        for lm_iter = 1:30
                            % LM update via SVD: delta = -V * (Sigma^2 + lambda^2*I)^{-1} * Sigma * U' * epsilon
                            d = sigma_vec ./ (sigma_vec.^2 + lambda^2);
                            delta_theta = -V_svd * (d .* (U_svd' * epsilon_real));
                            
                            % Update and project back to unit sphere
                            theta_try = theta + delta_theta;
                            theta_try = theta_try / norm(theta_try);  % Normalize to ||theta||=1
                            
                            % Evaluate new weighted cost
                            B_try = theta_try(1:nb+1);
                            A_try = theta_try(nb+2:end);
                            B_z_try = Z_b * B_try;
                            A_z_try = Z_a * A_try;
                            G_try = B_z_try ./ A_z_try;
                            epsilon_try = W_vec .* (G - G_try);  % Weighted error
                            cost_new = sum(abs(epsilon_try).^2);
                            
                            if cost_new < cost_old
                                % Step accepted
                                lambda = max(lambda / lambda_down, lambda_min);
                                theta = theta_try;
                                cost_old = cost_new;
                                step_accepted = true;
                                
                                if cost_new < best_cost
                                    best_cost = cost_new;
                                    best_theta = theta_try;
                                end
                                break;
                            else
                                % Step rejected
                                lambda = min(lambda * lambda_up, lambda_max);
                            end
                        end
                        
                        % Convergence checks (compare to cost BEFORE this iteration)
                        if step_accepted && abs(cost_new - cost_before_iter) / (cost_before_iter + eps) < 1e-10
                            n_lm_iter = iter;
                            break;
                        end
                        
                        if norm(delta_theta) / (norm(theta) + eps) < 1e-12
                            n_lm_iter = iter;
                            break;
                        end
                        
                        if lambda >= lambda_max
                            n_lm_iter = iter;
                            break;
                        end
                        
                        n_lm_iter = iter;
                    end
                    
                    % Use best solution found
                    theta = best_theta;
                    n_iter = n_lm_iter;
                    
                    % Extract B and A (full parameterization)
                    B = theta(1:nb+1).';
                    A = theta(nb+2:end).';
                    
                    % Normalize so a0 = 1 for output consistency
                    B = B / A(1);
                    A = A / A(1);
                    
                    fprintf('ML: init_cost=%.6e, final_cost=%.6e, improvement=%.2f%%, iters=%d\n', ...
                        cost_init, best_cost, 100*(1-best_cost/cost_init), n_lm_iter);
    end
    
    % Evaluate model
    z_all = exp(1j*2*pi*f/fs);
    G_model = (z_all.^(-(0:nb)) * B.') ./ (z_all.^(-(0:na)) * A.');
end

function theta = solve_ls(Phi, Y, w)
    W = diag(sqrt(w));
    Phi_w = [real(W*Phi); imag(W*Phi)];
    Y_w = [real(W*Y); imag(W*Y)];
    theta = Phi_w \ Y_w;
end
