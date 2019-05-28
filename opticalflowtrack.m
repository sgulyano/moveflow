function [ swc, flag ] = opticalflowtrack( swc, U, V )
%OPTICALFLOWTRACK tracking neuron using only optical flow

dx = interp2(U, swc(:,3), swc(:,4));
dy = interp2(V, swc(:,3), swc(:,4));

dx(isnan(dx)) = 0;
dy(isnan(dy)) = 0;

n = size(swc,1);
swc = swc + [zeros(n,2) dx dy zeros(n,3)];

s = size(U);
maxsize = ones(n,1)*s([2 1]);
count_out = sum(any(swc(:,3:4) < 1 | swc(:,3:4) > maxsize, 2));

% check if half of neuron is out-of-bound, then tracking is done
if count_out > n/2
    flag = true;
else
    flag = false;
end

