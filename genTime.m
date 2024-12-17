function [tList, dt] = genTime(tsteps, Tend)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if(nargin == 0)
    Tend = 8.0;
    tsteps = 1000000;
end
dt = Tend/tsteps;
fprintf('Time step: %f\n', dt);
% list of time steps
tList = 0:dt:Tend;
end