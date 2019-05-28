function [Oz, Sbest, E, Es] = estimate_plane_z(config, dx, Oz_prev)
%ESTIMATE_PLANE_Z esimate the depth of the neuron fitting plane

ang1 = cumsum(config{1});
ang2 = cumsum(config{2});

Ox = {-cos(ang1)*dx, cos(ang2)*dx};
Oz = {abs(sin(ang1)*dx), abs(sin(ang2)*dx)};
Sinit = {ones(size(Ox{1})), ones(size(Ox{2}))};

%% annealing to find the opimize zflip based on depth and curvature
Scurrent = Sinit;
[Ecurrent, Es] = ffd_depth_energy(Ox, Oz, Scurrent, Oz_prev);
fprintf('Init : E=%f Edep=%f Ecur=%f Eprev=%f\n', Ecurrent, Es(1), Es(2), Es(3));
Sbest = Scurrent; Ebest = Ecurrent;
maxiter = 100;
max_temp = 20;
temp_change = 0.95;
temp = max_temp;
for iter = 1:maxiter
    temp = temp * temp_change;
    Snew = Scurrent;
    % perturb
    gr = randi(2);
    i = randi(length(Ox{gr}));
    
    for dstate_i = [-1 1]
        Snew{gr}(i) = dstate_i;
        Enew = ffd_depth_energy(Ox, Oz, Snew, Oz_prev);
        
        if Enew <= Ecurrent || exp( (Ecurrent-Enew)/temp ) > rand()
            Scurrent = Snew;
            Ecurrent = Enew;
        end
        if Enew <= Ebest
            Sbest = Snew;
            Ebest = Enew;
        end
    end
end
[E, Es] = ffd_depth_energy(Ox, Oz, Sbest, Oz_prev);
fprintf('E=%f Edep=%f Ecur=%f Eprev=%f\n', E, Es(1), Es(2), Es(3));


end

function [E, Es] = ffd_depth_energy(Ox, Oz, S, Oz_prev)
Ox_all = [fliplr(cumsum(Ox{1})) 0 cumsum(Ox{2})];
Oz_all = [fliplr(cumsum( Oz{1}.*S{1}  )) 0 cumsum( Oz{2}.*S{2} )];

% Edep = sum(exp(Oz_all - min(Oz_all)));
Edep = sum((Oz_all - min(Oz_all)).^2);
Ecur = sum(diff(Ox_all,2).^2 + diff(Oz_all,2).^2);
Eprev = sum((Oz_all - Oz_prev).^2);

lines = [Ox_all(1:end-1)' Oz_all(1:end-1)' Ox_all(2:end)' Oz_all(2:end)'];
intrsect = false(size(lines,1));
% for i = size(lines,1)-11:size(lines,1)
for i = 1:size(lines,1)
    intrsect(i,:) = line_seg_intersect(lines(i,:), lines);
end

if any(any(tril(intrsect,-2) | triu(intrsect,2)))
    % self-intersection occurs
    E = inf;
    Es = [inf inf inf];
else
    E = Edep/2 + Ecur/4 + Eprev/4;
    Es = [Edep/2 Ecur/4 Eprev/4];
end
end

function out = line_seg_intersect(l1, l2)
% https://en.wikipedia.org/wiki/Intersection_(Euclidean_geometry)
tol = 0.5;
float_tol = 1e-10;

dx1 = l1(:,3) - l1(:,1);
dy1 = l1(:,4) - l1(:,2);
dx2 = l2(:,3) - l2(:,1);
dy2 = l2(:,4) - l2(:,2);
cx = l2(:,1) - l1(:,1);
cy = l2(:,2) - l1(:,2);

parallel_idx = abs(dx1*dy2 - dx2*dy1) < float_tol;
s = (cy.*dx2 - cx.*dy2) ./ (dx2.*dy1 - dx1.*dy2);
t = (dx1*cy - dy1*cx) ./ (dx2.*dy1 - dx1.*dy2);
% for non-parallel
out1 = isfinite(s) & isfinite(t) & s >= 0 & s <= 1 & t >= 0 & t <= 1;
% for parallel
p1 = (cx*dx1 + cy*dy1) / sqrt(dx1^2+dy1^2);
d1 = sqrt(cx.^2 + cy.^2 - p1.^2);
p1abs = p1 / sqrt(dx1^2+dy1^2);

cx2 = l2(:,3) - l1(:,1);
cy2 = l2(:,4) - l1(:,2);
p2 = (cx2*dx1 + cy2*dy1) / sqrt(dx1^2+dy1^2);
d2 = sqrt(cx2.^2 + cy2.^2 - p2.^2);
p2abs = p2 / sqrt(dx1^2+dy1^2);

out2 = parallel_idx & (d1 <= tol | d2 <= tol) & ((p1abs >= -float_tol & p1abs <= 1+float_tol) | (p2abs >= -float_tol & p2abs <= 1+float_tol));
out = out1 | out2;
end
