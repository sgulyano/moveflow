function [E, Ireg, Es] = ffd_energy( I_bs, I, Tx, yrange, config, constrain )
%FFD_ENERGY find FFD energy

%% image energy
% interpolate to transform I at time t to baseline coordinate
% [xx, yy] = meshgrid(Tx, yrange(1):yrange(2));
% keyboard;

Ireg = zeros(size(I_bs));
Tx_floor = floor(Tx);
Tx_next = Tx_floor+1;
Tx_alpha = Tx - Tx_floor;
Ty = yrange(1):yrange(2);
idxXIn = Tx_floor >= 1 & Tx_next <= size(I,2);
idxYIn = Ty >= 1 & Ty <= size(I,1);
I1 = I(Ty(idxYIn),Tx_floor(idxXIn));
I2 = I(Ty(idxYIn),Tx_next(idxXIn));
Ireg(idxYIn,idxXIn) = (ones(sum(idxYIn),1)*(1-Tx_alpha(idxXIn))).*I1 + (ones(sum(idxYIn),1)*Tx_alpha(idxXIn)).*I2;

% Ireg = interp2(double(I), xx, yy, 'linear', 0);

% compute NCC (normalized cross-correlation)
Iregnorm = Ireg - mean(Ireg(:));
Iregnorm = Iregnorm / sqrt(sum(sum(Iregnorm.^2)));
I_bsnorm = I_bs - mean(I_bs(:));
I_bsnorm = I_bsnorm / sqrt(sum(sum(I_bsnorm.^2)));

Eimg = 1 - max(sum(sum(Iregnorm .* I_bsnorm)), 0);

%% smooth energy
dtheta = cellfun(@(x)abs(diff(x,[],2)), config, 'UniformOutput', false);
dtheta = [fliplr(dtheta{1}) config{1}(1)+config{1}(2) dtheta{2}];
dtheta(dtheta > pi) = 2*pi - dtheta(dtheta > pi);

% dtheta = config(1:end-1) - config(2:end);
% % at soma, flip dtheta
% dtheta(ctrlpnt_dist(1)) = config(ctrlpnt_dist(1)) + config(ctrlpnt_dist(1)+1);
% dtheta = abs(dtheta);
% dtheta(dtheta > pi) = 2*pi - dtheta(dtheta > pi);

Esmooth = sum(dtheta.^2);

%% shape energy that keep dendrites flat
theta = [config{1} config{2}];
Eshape = sum(theta.^2);

%% constraints energy at dendrite tips and soma
Econ = sum((Tx(constrain(:,1)) - constrain(:,2)').^2);

%% total energy
E = 10*Eimg + Esmooth + 20*Eshape + Econ;
Es = {10*Eimg, Esmooth, Eshape, Econ};

end

