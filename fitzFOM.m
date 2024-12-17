function Y = fitzFOM(N, tsteps, Tend, prefix, storeI)

% use empty prefix, if no prefix is given
if(nargin < 4)
    prefix = '';
end

% use default tstaps/1000 for taking snapshots
if(nargin < 5)
    mys = 1:tsteps;
    storeI = mys(mod(mys, tsteps/1000) == 0);
end

% system parameters and operators
[h, e, A, c] = genOperators(N);
% time step
[tList, dt] = genTime(tsteps,Tend);

% Initialize vectors with initial conditions y(0)=0.
y_previous = zeros(2*N, 1);
Y = [];
for i=2:tsteps+1
    % compute vector g + c + f
    g_plus_c_plus_f = c;
    % f affects only the upper 1:N-entries, since 
    % f = [nonlin(y_previous(1:N)); zeros(N, 1)];
    g_plus_c_plus_f(1:N) = g_plus_c_plus_f(1:N)+(1/e)*nonlin(y_previous(1:N));
    % g affecst only the first entry, i.e. boundary initial condition
    g_plus_c_plus_f(1) = g_plus_c_plus_f(1) + e/h*bc(tList(i));
    
    % set up right hand side
    rhs = A*y_previous + g_plus_c_plus_f;
    % compute next y
    y_new = dt*rhs + y_previous;
    %reset
    y_previous = y_new;
    
    % check, if we take a snapshot
    if(find(storeI == i, 1))
        disp(tList(i));
        Y(:, end + 1) = y_new;
    end
end
% store the last snapshot
Y(:, end + 1) = y_new;

eval(['save snapshots_fitz/', prefix, 'snapshots_N', num2str(N), '_tsteps', num2str(tsteps), '_Tend', num2str(Tend), '.mat Y tsteps storeI']);

end


