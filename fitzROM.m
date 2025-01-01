function Y_reconstructed = fitzROM(N, tsteps, Tend, prefix, storeI)

% Use empty prefix, if no prefix is given
if(nargin < 4)
    prefix = '';
end

if(nargin < 5)
    % Save every 1000th time step
    mys = 1:tsteps;
    storeI = mys(mod(mys, tsteps/1000) == 0);
end


filePath = sprintf('snapshots_fitz/%ssnapshots_N%d_tsteps%d_Tend%d.mat', prefix, N, tsteps, Tend);

snapshots = load(filePath); % Load snapshots
threshold = 99.99; % Threshold for RIC; change as desired.

Y = snapshots.Y; % Size 2N x M
N = size(Y, 1) / 2; % Number of rows in V and W
% Extract V and W from snapshots
V = Y(1:1024, :);   % First 1024 rows
W = Y(1025:end, :); % Last 1024 rows

% Parameter and operator setup
[h, e, A, c] = genOperators(N);
[tList, dt] = genTime(tsteps,Tend);

% ------------------------- SVD & MODES -------------------------
[U_V, S_V, ~] = svd(V, 'econ');
[U_W, S_W, ~] = svd(W, 'econ');

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
hold off;

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


prompt = "Number of modes to retain. Empty input will use RIC. Modes: ";
modes = input(prompt);
tic; % As we are no longer inflating task time by waiting for input

if(not(isempty(modes)) & isnumeric(modes))    
    p = modes;
end



k_V = PODModes(diag(S_V), threshold);
k_W = PODModes(diag(S_W), threshold);
if(isempty(modes) | not(isnumeric(modes)))
    fprintf('Number of modes retained for V: %d\n', k_V);
    fprintf('Number of modes retained for W: %d\n', k_W);
    p = max(k_V, k_W); % Ensures greatest accuracy
end

fprintf('Using modes p = %f\n', p)


% ------------------- REDUCED ORDER MODEL -------------------
% Get the reduced bases
U_V = U_V(:,1:p);
U_W = U_W(:,1:p);
U = [U_V, zeros(N, p); zeros(N, p), U_W];  % Size: 2N x 2p


% ------------------- POD & EXPLICIT EULER ---y----------------
% Project the constant full-order operators onto the reduced basis U
A_r = U' * A * U;
c_r = U' * c;

% Define e1, the first basis vector
e_1 = zeros(N, 1);
e_1(1) = 1;

% Initial conditions
y_previous = zeros(2*p, 1); 

% Initialize g, F; after p rows (i.e. when w), they are 0.
g_r = zeros(2*p, 1);
F_r = zeros(2*p, 1);

% Time-stepping loop
Y_snapshots = [];  % Initialize snapshot matrix
for i = 2:tsteps+1

    % Compute g, F; for the first p rows (i.e. when v), they are not necessarily 0.
    g_r(1:p) = e/h * bc(tList(i)) * (U_V' * e_1);
    F_r(1:p) = (1/e) * U_V' * nonlin(U_V * y_previous(1:p));

    % Compute g + F + c
    gFc_sum = g_r + F_r + c_r;

    % Compute the next reduced-order solution
    y_new = y_previous + dt * (A_r * y_previous + gFc_sum);
    y_previous = y_new;

    % Check if we take a snapshot
    if any(storeI == i)
        disp(tList(i));
        Y_snapshots(:, end + 1) = y_new;
    end
end
% Store the last snapshot
Y_snapshots(:, end + 1) = y_new;

% Reconstruct the reduced-order solution
Y_reconstructed = U * Y_snapshots;
V_reconstructed = Y(1:p,:);
W_reconstructed = Y(p+1:end,:);

err_Y = Y - Y_reconstructed;
fro_err_Y = norm(err_Y, 'fro');
fro_norm_Y = norm(Y, 'fro');
rel_err_Y = fro_err_Y / fro_norm_Y;

% Display the results
fprintf('Frobenius Norm of the Error: %.4e\n', fro_err_Y);
fprintf('Relative Error: %.4e\n', rel_err_Y);
toc;

end
