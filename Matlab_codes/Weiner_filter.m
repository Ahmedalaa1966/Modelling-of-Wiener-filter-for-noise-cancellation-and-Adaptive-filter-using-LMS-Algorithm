clc;
clear;
close all;

%% Parameters
N = 1000;
w0 = 0.05*pi;
phi = 0;
orders = [6 12];

% Generate desired signal
n = 0:N-1;
d = sin(w0*n + phi);

%% Generate white noise
g = randn(1,N);
g = g - mean(g);      % zero mean
g = g / std(g);       % unit standard deviation


v1 = zeros(1,N);
for k = 2:N
    v1(k) = 0.8*v1(k-1) + g(k);       % starting the counter form 2 to avoid negative index
end

v2 = zeros(1,N);
for k = 2:N                            % starting the counter form 2 to avoid negative index
    v2(k) = -0.6*v2(k-1) + g(k);
end

x = d + v1;                            % Corrupted signal

% Wiener filter design
for idx = 1:length(orders)             % for loop that loops twice one for order 12 and one for order 6

    M = orders(idx);
    % Estimate autocorrelation of v2 (rv2)
    rv2 = zeros(M+1,1);

    for k = 0:M
        sum_val = 0;
        for n1 = k+1:N
            sum_val = sum_val + v2(n1)*v2(n1-k);      % no complex conjuct as the function is real
        end
        rv2(k+1) = sum_val/N;
    end

    % Form autocorrelation matrix
    Rv2 = toeplitz(rv2);

    % Estimate cross-correlation
    rxv2 = zeros(M+1,1);

    for k = 0:M

        sum_val = 0;

        for n1 = k+1:N
            sum_val = sum_val + x(n1)*v2(n1-k);
        end

        rxv2(k+1) = sum_val/N;
    end

    %% Wiener solution
    w = inv(Rv2)*rxv2;

    %% Filter output
    v1_hat = filter(w,1,v2);

    %% Error signal
    e = x - v1_hat;

    %% Plot
    figure;

    subplot(3,1,1);
    plot(n,x);
    title(['Corrupted Signal x(n), Order = ',num2str(M)]);

    subplot(3,1,2);
    plot(n,e);
    title('Recovered Signal e(n)');

    subplot(3,1,3);
    plot(n,d);
    title('Original Signal d(n)');

    %% MSE
    mse = mean((d-e).^2);

    disp(['Filter Order = ',num2str(M)]);
    disp(['MSE = ',num2str(mse)]);
end