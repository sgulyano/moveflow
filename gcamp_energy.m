function [E, Ireg, Es] = gcamp_energy( I_bs, I, Tx, dY )
%GCAMP_ENERGY Summary of this function goes here
%   Detailed explanation goes here


%% image energy
% interpolate to transform I at time t to baseline coordinate
% [xx, yy] = meshgrid(Tx, yrange(1):yrange(2));
% keyboard;

Ireg = zeros(size(I_bs));
Tx_floor = floor(Tx);
Tx_next = Tx_floor+1;
Tx_alpha = Tx - Tx_floor;
Ty = (1:size(I_bs,1)) + round(dY);
idxXIn = Tx_floor >= 1 & Tx_next <= size(I,2);
idxYIn = Ty >= 1 & Ty <= size(I,1);
I1 = I(Ty(idxYIn),Tx_floor(idxXIn));
I2 = I(Ty(idxYIn),Tx_next(idxXIn));
Ireg(idxYIn,idxXIn) = (ones(sum(idxYIn),1)*(1-Tx_alpha(idxXIn))).*I1 + (ones(sum(idxYIn),1)*Tx_alpha(idxXIn)).*I2;

% % Ireg = interp2(double(I), xx, yy, 'linear', 0);
% 
% % compute NCC (normalized cross-correlation)
% Iregnorm = Ireg - mean(Ireg(:));
% Iregnorm = Iregnorm / sqrt(sum(sum(Iregnorm.^2)));
% I_bsnorm = I_bs - mean(I_bs(:));
% I_bsnorm = I_bsnorm / sqrt(sum(sum(I_bsnorm.^2)));
% 
% Eimg = 1 - max(sum(sum(Iregnorm .* I_bsnorm)), 0);

Eimg = sum(sum((Ireg - I_bs).^2));

%% smooth energy
dTx = diff(Tx);
Esmooth = sum((dTx - 1).^2);

%% total energy
E = Eimg/20 + Esmooth;
Es = {Eimg/20, Esmooth};


end

