%% Task 2: Adaptive Linear Prediction Using the LMS Algorithm
% AR(2) process: x(n) = 1.2728*x(n-1) - 0.81*x(n-2) + v(n)
% Optimum predictor coefficients: w1 = 1.2728, w2 = -0.81
% sigma_v^2 = 1

clear; close all; clc;

%% =========================================================
%  PART 1: Generate the AR(2) signal x(n)
%% =========================================================
N = 1000;
v = randn(N, 1);   % unit-variance white noise

x = zeros(N, 1);
for n = 3:N
    x(n) = 1.2728*x(n-1) - 0.81*x(n-2) + v(n);
end

%% =========================================================
%  PART 2: LMS Adaptive Linear Predictor
%  x_hat(n) = w_n(1)*x(n-1) + w_n(2)*x(n-2)
%  Update: w_{n+1}(k) = w_n(k) + mu * e(n) * x(n-k)
%% =========================================================
mu_values = [0.02, 0.04];
colors    = {'b', 'r'};

figure('Name','Convergence of Filter Coefficients','NumberTitle','off');

for idx = 1:length(mu_values)
    mu = mu_values(idx);

    w      = zeros(2, 1);        % initial weights = 0
    W_hist = zeros(N, 2);
    E_sq   = zeros(N, 1);

    for n = 3:N
        x_vec = [x(n-1); x(n-2)];
        x_hat = w' * x_vec;
        e     = x(n) - x_hat;
        w     = w + mu * e * x_vec;  % LMS update

        W_hist(n, :) = w';
        E_sq(n)      = e^2;
    end

    % Plot: Coefficient Convergence
    subplot(2, 2, (idx-1)*2 + 1);
    plot(W_hist(:,1), 'b', 'LineWidth', 1.2); hold on;
    plot(W_hist(:,2), 'r', 'LineWidth', 1.2);
    yline(1.2728, 'b--', 'LineWidth', 1);
    yline(-0.81,  'r--', 'LineWidth', 1);
    xlabel('Iteration'); ylabel('Coefficients');
    title(sprintf('Coefficients Convergence  (\\mu = %.2f)', mu));
    legend('w_n(1)', 'w_n(2)', 'w_1^{opt}=1.2728', 'w_2^{opt}=-0.81', ...
           'Location','east');
    grid on; xlim([1 N]);

    % Plot: Squared Prediction Error
    subplot(2, 2, (idx-1)*2 + 2);
    plot(E_sq, colors{idx}, 'LineWidth', 0.8);
    xlabel('Iteration'); ylabel('Squared Error');
    title(sprintf('Squared Prediction Error  (\\mu = %.2f)', mu));
    grid on; xlim([1 N]);
    if mu == 0.04
        ylim([0 20]);
    end
end
sgtitle('LMS Adaptive Linear Predictor — Convergence');

%% =========================================================
%  PART 3: PSD via Z-transform, then residuez for r_x(k)
%
%  H(z) = 1 / (1 - 1.2728z^{-1} + 0.81z^{-2})
%  P_x(z) = sigma_v^2 * H(z) * H(1/z)
%
%  H(1/z) = 1/(1-1.2728z+0.81z^2) = z^{-2}/(0.81-1.2728z^{-1}+z^{-2})
%
%  So P_x(z) = z^{-2} / [(1-1.2728z^{-1}+0.81z^{-2})(0.81-1.2728z^{-1}+z^{-2})]
%
%  Use residuez -> partial fractions -> inverse Z-transform gives r_x(k)
%
%  KEY RULE: poles with |p|<1 are CAUSAL -> r_x(k>=0) = sum R_i * p_i^k
%            poles with |p|>1 are ANTI-CAUSAL -> contribute only for k<0
%% =========================================================

%% 3a) Define PSD as ratio of polynomials in z^{-1}
num_psd = [0, 0, 1];                           % z^{-2} numerator
den_psd = conv([1, -1.2728, 0.81], ...         % den of H(z)
               [0.81, -1.2728, 1]);             % den of H(1/z) in z^{-1}

%% 3b) Partial fraction expansion
[R_pf, p, C] = residuez(num_psd, den_psd);

fprintf('\n--- Partial Fraction Expansion (residuez) ---\n');
fprintf('Poles:\n'); disp(p);
fprintf('Residues:\n'); disp(R_pf);
fprintf('|Poles|:\n'); disp(abs(p));

%% 3c) Extract CAUSAL poles only (|p| < 1) for r_x(k >= 0)
%
%  For the two-sided PSD:
%    Causal poles   (|p|<1)  -> r_x(k) for k >= 0  via  R_i * p_i^k
%    Anti-causal poles (|p|>1) -> r_x(k) for k <  0  (not needed here)
%
%  Therefore:  r_x(k) = real( sum_i  R_i * p_i^k )   for k >= 0
%              where the sum is over CAUSAL poles ONLY

causal_mask = abs(p) < (1 - 1e-9);
p_c = p(causal_mask);
R_c = R_pf(causal_mask);

fprintf('\n--- Causal poles (|p| < 1) ---\n');
disp(p_c);

% Compute r_x(0), r_x(1), r_x(2)
r0 = real(sum(R_c));                % p_i^0 = 1 for all
r1 = real(sum(R_c .* p_c.^1));
r2 = real(sum(R_c .* p_c.^2));

fprintf('\n--- Autocorrelation (analytical via residuez, causal poles only) ---\n');
fprintf('r_x(0) = %.4f   (lecture answer: 5.7524)\n', r0);
fprintf('r_x(1) = %.4f   (lecture answer: 4.0451)\n', r1);
fprintf('r_x(2) = %.4f\n', r2);

%% 3d) Plot the PSD (analytical)
omega = linspace(-pi, pi, 1024);
z_eval = exp(1j*omega);
H_z    = 1 ./ (1 - 1.2728*z_eval.^(-1) + 0.81*z_eval.^(-2));
Pxx_analytical = abs(H_z).^2;      % sigma_v^2 = 1

figure('Name','Power Spectral Density (Analytical)','NumberTitle','off');
plot(omega/pi, 10*log10(Pxx_analytical), 'b', 'LineWidth', 1.5);
xlabel('Normalised Frequency (\times\pi rad/sample)');
ylabel('PSD (dB)');
title('Power Spectral Density of x(n) — Analytical');
grid on;

%% 3e) Autocorrelation matrix R (2x2 Toeplitz)
R_mat = toeplitz([r0, r1]);

fprintf('\n--- Autocorrelation Matrix R ---\n');
fprintf('R = [%.4f  %.4f]\n', R_mat(1,1), R_mat(1,2));
fprintf('    [%.4f  %.4f]\n', R_mat(2,1), R_mat(2,2));
fprintf('(Lecture answer: R = [5.7524  4.0451; 4.0451  5.7524])\n');

%% 3f) Eigenvalues and step-size upper bound
eigenvalues = eig(R_mat);
lambda_max  = max(real(eigenvalues));
lambda_min  = min(real(eigenvalues));

% LMS convergence bound: 0 < mu < 2 / lambda_max
mu_max = 2 / lambda_max;

fprintf('\n--- Eigenvalues of R ---\n');
fprintf('lambda_max = %.4f   (lecture answer: 9.7975)\n', lambda_max);
fprintf('lambda_min = %.4f\n', lambda_min);
fprintf('\n--- Step-size Upper Bound:  0 < mu < 2/lambda_max ---\n');
fprintf('mu_max = 2 / %.4f = %.4f   (lecture answer: 0.2041)\n', lambda_max, mu_max);
fprintf('\nNote: mu should be ~1 order of magnitude below mu_max.\n');
fprintf('      mu=0.02 is a good choice (lecture recommendation).\n');
fprintf('      Both 0.02 and 0.04 satisfy the bound.\n');

%% =========================================================
%  PART 4: Final Summary
%% =========================================================
fprintf('\n============= SUMMARY =============\n');
fprintf('Optimum coefficients : w1 = 1.2728,  w2 = -0.81\n');
fprintf('Minimum MSE          : xi_min = sigma_v^2 = 1\n');
fprintf('r_x(0) = %.4f,  r_x(1) = %.4f\n', r0, r1);
fprintf('R = [%.4f  %.4f;\n     %.4f  %.4f]\n', R_mat(1,1),R_mat(1,2),R_mat(2,1),R_mat(2,2));
fprintf('lambda_max = %.4f,  lambda_min = %.4f\n', lambda_max, lambda_min);
fprintf('Step-size bound: 0 < mu < %.4f\n', mu_max);