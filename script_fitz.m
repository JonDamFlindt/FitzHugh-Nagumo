clear;
close all;

% Before you run this script, execute the function fitzFOM
% tic; fitzFOM(1024, 1000000, 8.0, 'Test_'); toc;
% This may take a while.
% You might want to record the computation time.
%
% Afterwards, edit the parameters below in FitzROM to match.
% tic; fitzROM(1024, 1000000, 8.0, 'Test_'); toc;
% This shouldn't take that long as long as the amount of modes p is sensible.
% You might want to record the computation time.
%
%

% if your files/folders have different names, modify accordingly
Snaps = load('snapshots_fitz/Test_snapshots_N1024_tsteps1000000_Tend8.mat');
Y = Snaps.Y;

% use the same parameters as the ones loaded into Snaps
Y_reconstructed = fitzROM(1024, 1000000, 8.0, 'Test_');


t = linspace(0, 8.0, size(Y, 2)); % Linear time space for the x-axis.
V = Snaps.Y(1:1024,:);
W = Snaps.Y(1025:end,:);

% ------------------- REDUCED ORDER MODEL -------------------
V_reconstructed = Y_reconstructed(1:1024, :);
W_reconstructed = Y_reconstructed(1025:end, :);

% Generate time vector for snapshots
t_reconstructed = linspace(0, 8.0, size(Y_reconstructed,2));

figure;
plot(t, V(1, :), 'r-');
hold on;
plot(t, W(1, :), 'b-');
xlabel('Time');
ylabel('Voltage');
legend('V', 'W');
title('Full-order model solutions');

% Plot both solutions
figure;
plot(t, V(1, :), 'r-');
hold on;
plot(t, W(1, :), 'r-');
plot(t_reconstructed, V_reconstructed(1, :), 'b--');
plot(t_reconstructed, W_reconstructed(1, :), 'b--');
xlabel('Time');
ylabel('Voltage');
legend('V', 'W', 'V reconstructed', 'W reconstructed');
title('Reduced-order model solutions');
hold off;

