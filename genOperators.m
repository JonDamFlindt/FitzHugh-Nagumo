function [h, e, A, c] = genOperators(N)

L = 1; % length
e = 0.015;
b = 0.5;
gamma = 2;
c = 0.05;
h = L/(N + 2);

Dxx = toeplitz([-2 1 zeros(1, N - 2)]);
Dxx(1, 1) = -1;
Dxx(end, end) = -1;
Dxx = (1/h^2)*Dxx;
A = [e*Dxx, (-1/e)*eye(N); b*eye(N), -gamma*eye(N)];
c = c*ones(2*N, 1);
c(1:N) = (1/e)*c(1:N);
return;
end
