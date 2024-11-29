function [E, A, c] = genOperators(N)

L = 1; % length
e = 0.015;
b = 0.5;
gamma = 2;
c = 0.05;
h = L/(N + 2);

E = ones(2*N,1);
E(1:N) = e*E(1:N);
% Caution: this E does not have the factor (1/dt) !!
K = toeplitz([-2 1 zeros(1, N - 2)]);
K(1, 1) = -1;
K(end, end) = -1;
K = -K;
A = [-e^2/h^2*K, -eye(N); b*eye(N), -gamma*eye(N)];
c = c*ones(2*N, 1);
return;
end