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


% ------------------ FULL ORDER MODEL ------------------
figure;
hold on;
plot(t, V(1,:));
plot(t, W(1,:));
xlabel('Time Step');
ylabel('Voltage');
title('First Row of Voltage V and Recovery Voltage W');
legend('V', 'W')
hold off;


% ------------------------- SVD -------------------------
[U_V, S_V, ~] = svd(V, 'econ');
[U_W, S_W, ~] = svd(W, 'econ');


k_V = PODModes(diag(S_V));
fprintf('Number of modes retained for V: %d\n', k_V);
k_W = PODModes(diag(S_W));
fprintf('Number of modes retained for W: %d\n', k_W);
p = min(k_V, k_W); % Should this be max(k_V, k_W)? Neither/independent?



% ------------------- REDUCED ORDER MODEL -------------------
% Get the reduced bases
U_V = U_V(:,1:p);
U_W = U_W(:,1:p);
U = [U_V, zeros(N, p); zeros(N, p), U_W];  % Size: 2N x (k_V + k_W)

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


% ------------------- POD & EXPLICIT EULER -------------------
[E, A, c] = genOperators(N);

% Project the constant full-order operators onto the reduced basis U
A_r = U' * A * U;
c_r = U' * c;

% Parameter setup
Tend = 8.0;
tsteps = 1000000;
dt = Tend/tsteps;
fprintf('Time step: %f\n', dt);
% list of time steps
tList = 0:dt:Tend;
L = 1; % length
% system parameters:
e = 0.015;
h = L/(N + 2);

% Save every 1000th time step
times = 1:tsteps;
storeInterval = times(mod(times, tsteps/1000) == 0);

% 90% sure that the error is below this line
%---------------------------------------- 

% Define e1
% What in the world is e1? I guess it's the first basis vector?
e_1 = zeros(N, 1);
e_1(1) = 1;


% Define the initial conditions, y(0)=0
y_previous = zeros(2*p, 1); % Correct, for all spatial points.


% Time-stepping loop
Y = [];  % Initialize snapshot matrix
for k = 2:tsteps+1
    t_k = tList(k);

    % Initialize g, F; after p rows (i.e. when w), they are 0.
    g_r = zeros(2*p, 1);
    F_r = zeros(2*p, 1);

    % Compute g, F; for the first p rows (i.e. when v), they are not necessarily 0.
    g_r(1:p) = e/h * bc(t_k) * U_V' * e_1;
    F_r(1:p) = (1/e) * U_V' * nonlin(U_V * y_previous(1:p));

    % Compute g + F + c
    gFc_sum = g_r + F_r + c_r;

    % Compute the next reduced-order solution
    y_new = y_previous + dt * (A_r * y_previous + gFc_sum);
    y_previous = y_new;

    % Check if we take a snapshot
    if any(storeInterval == k)
        disp(t_k);
        Y(:, end + 1) = y_new;
    end
end

% Reconstruct the reduced-order solution
Y_full = U * Y;
V_reconstructed = Y_full(1:N, :);
W_reconstructed = Y_full(N+1:end, :);

% Generate time vector for snapshots
t_snapshots = tList(storeInterval);

% Plot the solutions
figure;
plot(t_snapshots, V_reconstructed(1, :));
hold on;
plot(t_snapshots, W_reconstructed(1, :));
xlabel('Time');
ylabel('Solution (Voltage)');
legend('V', 'W');
title('Reduced-order model solutions');
hold off;

% Plot the solutions
figure;
plot(W_reconstructed(1,:), V_reconstructed(1,:), 'r--');
hold on;
%plot(tList, W_reconstructed(1, :));
xlabel('W');
ylabel('V');
title('Phase plot of V and W');
