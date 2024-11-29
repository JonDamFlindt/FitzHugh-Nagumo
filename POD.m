clear; close all;
% Load snapshots
snapshots = load('snapshots_fitz/Test_snapshots_N1024_tsteps1000000_Tend8.mat'); % Adjust file path as needed
Y = snapshots.Y; % Size 2N x M

Tend = 8.0;
num_snapshots = size(Y, 2);
t = linspace(0, Tend, num_snapshots);

N = size(Y, 1) / 2; % Number of rows in V and W
% Extract V and W from snapshots
V = Y(1:1024, :);   % First 1024 rows
W = Y(1025:end, :); % Last 1024 rows


figure;
hold on;
plot(t, V(1,:));
plot(t, W(1,:));
xlabel('Time Step');
ylabel('Voltage');
title('First Row of Voltage V and Recovery Voltage W');
legend('V', 'W')
hold off;


% SVD
[U_V, S_V, ~] = svd(V, 'econ');
[U_W, S_W, ~] = svd(W, 'econ');


function k = GetPODModes(singular_values)
    % Change as desired, personally I prefer 100 - 1e-6
    threshold = 100 - 1.0e-6;
    
    % Calculate the Relative Info Content
    info_sum = sum(singular_values);
    info_cumsum = cumsum(singular_values);
    Relative_Information_Content = (info_cumsum / info_sum) * 100;

    % Choose the first k that satisfies the threshold
    k = find(Relative_Information_Content >= threshold, 1);
end

k_V = GetPODModes(diag(S_V));
fprintf('Number of modes retained for V: %d\n', k_V);
k_W = GetPODModes(diag(S_W));
fprintf('Number of modes retained for W: %d\n', k_W);



% Get the reduced bases
U_V = U_V(:,1:k_V);
U_W = U_W(:,1:k_W);
U = [U_V, zeros(N, k_W); zeros(N, k_V), U_W];  % Size: 2N x (k_V + k_W)

% Plot singular values and RIC
figure;
subplot(2,1,1);
semilogy(diag(S_V), 'b', 'LineWidth', 2);
hold on;
semilogy(diag(S_W), 'r', 'LineWidth', 2);
xlabel('Number of columns');
ylabel('Singular values');
legend('V', 'W');
title('Singular values vs. Number of columns');

% Calculate RIC
RIC_V = cumsum(diag(S_V)) / sum(diag(S_V)) * 100;
RIC_W = cumsum(diag(S_W)) / sum(diag(S_W)) * 100;

subplot(2,1,2);
semilogx(RIC_V, 'b', 'LineWidth', 2);
hold on;
semilogx(RIC_W, 'r', 'LineWidth', 2);
xlabel('Number of columns');
ylabel('Relative Information Content (%)');
legend('V', 'W');
title('RIC vs. Number of columns');
hold off;

% Generate the operators
[E, A, c] = genOperators(N);

% Project the constant full-order operators onto the reduced basis U
A_r = U' * A * U;
c_r = U' * c;

% Parameter setup
Tend = 8.0;
tsteps = 125000; % Number of time steps, kept crash at 1000000
dt = Tend/tsteps;
fprintf('Time step: %f\n', dt);
% list of time steps
tList = 0:dt:Tend;
L = 1; % length
% system parameters:
e = 0.015;
h = L/(N + 2);

% 90% sure that the error is below this line
%---------------------------------------- 

% Define e1
% What in the world is e1? I guess it's the first basis vector?
e_1 = zeros(N, 1);
e_1(1) = 1;


% Assumption regarding the initial conditions
a0 = zeros(k_V + k_W, 1);  % Initial conditions
y = zeros(k_V + k_W, length(tList));
y(:, 1) = a0;
y(1:k_V) = -bc(0) * U_V' * e_1;


% Time-stepping loop
for k = 1:tsteps
    t_k = tList(k);


    % Compute g & F
    g_r = zeros(k_V + k_W, 1);
    F_r = zeros(k_V + k_W, 1);

    g_r(1:k_V) = g_r(1:k_V) + e/h * bc(t_k) * U_V' * e_1;
    F_r(1:k_V) = 1/e * U_V' * nonlin(U_V * y(1:k_V, k));

    % Compute g + F + c
    gFc_sum = g_r + F_r + c_r;

    % Compute the next reduced-order solution
    y(:, k+1) = y(:, k) + dt * (A_r * y(:, k) + gFc_sum);

    
end

% Reconstruct the reduced-order solution
Y_full = U * y;
V_reconstructed = Y_full(1:N, :);
W_reconstructed = Y_full(N+1:end, :);

% Plot the solutions
figure;
plot(tList, V_reconstructed(1, :));
hold on;
plot(tList, W_reconstructed(1, :));
xlabel('Time');
ylabel('Solution (Voltage)');
legend('V', 'W');
title('Reduced-order model solutions');
