function [E, Ireg, Es] = fit_plane_energy( I0, I1, S0, S1, xrange, yrange, dx, constrain )
%FFD_ENERGY find FFD energy

%% image energy
% interpolate to transform I at time t to baseline coordinate
% [xx, yy] = meshgrid(Tx, yrange(1):yrange(2));
% keyboard;

Ireg = zeros(size(I0));

Tx = interp1(S0, S1, xrange(1):xrange(2));
Tx_floor = floor(Tx);
Tx_next = Tx_floor+1;
Tx_alpha = Tx - Tx_floor;
Ty = yrange(1):yrange(2);
idxXIn = Tx_floor >= 1 & Tx_next <= size(I1,2);
idxYIn = Ty >= 1 & Ty <= size(I1,1);

I_floor = I1(Ty(idxYIn),Tx_floor(idxXIn));
I_next = I1(Ty(idxYIn),Tx_next(idxXIn));
Ireg(idxYIn,idxXIn) = (ones(sum(idxYIn),1)*(1-Tx_alpha(idxXIn))).*I_floor + (ones(sum(idxYIn),1)*Tx_alpha(idxXIn)).*I_next;
% Ireg = interp2(double(I), xx, yy, 'linear', 0);

% compute NCC (normalized cross-correlation)
Iregnorm = Ireg - mean(Ireg(:));
Iregnorm = Iregnorm / sqrt(sum(sum(Iregnorm.^2)));
I_bsnorm = I0 - mean(I0(:));
I_bsnorm = I_bsnorm / sqrt(sum(sum(I_bsnorm.^2)));

Eimg = 1 - max(sum(sum(Iregnorm .* I_bsnorm)), 0);

%% smooth energy
Esmooth = sum(diff(S1,2));

%% shape energy that keep dendrites flat
Eshape = sum(max(dx - diff(S1,1), 0)) + sum(exp(100*max(diff(S1,1) - dx, 0))-1);

%% constraints energy at dendrite tips and soma
anchor_pnt = interp1(S0, S1, constrain(:,1));
Econ = sum((anchor_pnt - constrain(:,2)).^2);

%% total energy
E = 10*Eimg + Esmooth + 20*Eshape + Econ;
Es = {10*Eimg, Esmooth, Eshape, Econ};
% keyboard;
end

