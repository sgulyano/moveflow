function [ newswc ] = adjustZ( swct, Vt, swcp, Vp, opt )
%ADJUSTZ adjust Z-coordinate of swct

if nargin < 5;  opt = struct();  end;
if ~isfield(opt,'alpha_d');       opt.alpha_d 	= 1/5;      end;
if ~isfield(opt,'alpha_t');     opt.alpha_t = 1/10;     end;
if ~isfield(opt,'reg');         opt.reg 	= 1/2;      end;
if ~isfield(opt,'r_oof');       opt.r_oof 	= 3;        end;
if ~isfield(opt,'DEBUG');       opt.DEBUG   = false;    end;


%% adjust Z using MRF over tree structure
% Generate Adjacency Matrix
MAXDEP = size(Vt,3);
states = 1:MAXDEP;
numLb = length(states);
numVar = size(swct,1);

% neuron graph
nb = swct(:,[1 7]);
nb = nb(all(nb>0,2),:);
T = sparse(nb(:,1), nb(:,2), 1, numVar, numVar);
T = T + T';

pnt = round(swct(:,3:4));


edgeStruct = UGM_makeEdgeStruct(T,numLb);

Vadj = Vt;
for ii = 1:size(Vt,3)
    Vadj(:,:,ii) = imadjust(Vt(:,:,ii), [.02, .4], [0 1]);
end
Vadj = im2double(Vadj);

%% find vesselness
% opts.spacing = [1, 1, 2.12 / 0.624];
% opts.responsetype = 4;
% Vpad = padarray(Vt, opt.r_oof+5*ones(1,3), 'both');
% Ivessel = oof3response(double(Vpad), 1:opt.r_oof, opts);
% % Options.BlackWhite = false;
% % Options.FrangiScaleRange = [1 2];
% % Options.FrangiScaleRatio = 0.5;
% % Options.verbose = false;
% % Options.FrangiAlpha = .1;
% % Options.FrangiBeta  = .9;
% % Options.FrangiC = 10000;
% % Ivessel = FrangiFilter3D(double(Vadj),Options);
% figure(2); imshow3D(Ivessel);

%% get linear index of nodes
p_x = pnt(:,1)*ones(1,numLb);
p_y = pnt(:,2)*ones(1,numLb);
p_z = ones(numVar,1)*states;
p_ind = sub2ind(size(Vt), p_y, p_x, p_z);

%% add unaries 
% Image feature: vesselness & function of eigenvalues of Hessian matrix
Eimg = Vadj(p_ind);%Ivessel(p_ind);

dz = bsxfun(@minus, swcp(:,5), states);
Emodel = 1.5*exp(-abs(dz) * opt.alpha_d);

nodePot = Eimg + Emodel;

%% add pairwise potential
edgePot = ones(numLb,numLb,edgeStruct.nEdges);
if opt.DEBUG
    Esmooths = zeros(numLb,numLb,edgeStruct.nEdges);
end
% Transition smoothness
dz = bsxfun(@minus, states, states');
dsmooth = abs(dz);
Esmooth = exp(-dsmooth * opt.alpha_t);
Esmooth = bsxfun( @rdivide, Esmooth, 2*mean(Esmooth(:)) );

for e = 1:edgeStruct.nEdges
    i = edgeStruct.edgeEnds(e,1);
    j = edgeStruct.edgeEnds(e,2);
    if T(i,j) == 1
        edgePot(:,:,e) = opt.reg * Esmooth;
        if opt.DEBUG
            Esmooths(:,:,e) = opt.reg * Esmooth;
        end
    else
        error('unknown edge type');
    end
end

%% inference
[~, optimalDecoding] = fusion_moves(-nodePot, -edgePot, edgeStruct.edgeEnds, 20);

nidx = sub2ind(size(nodePot),1:numVar,optimalDecoding');
Enode = nodePot(nidx);
eidx = sub2ind(size(edgePot),...
        optimalDecoding(edgeStruct.edgeEnds(:,1))',...
        optimalDecoding(edgeStruct.edgeEnds(:,2))',...
        1:edgeStruct.nEdges);
Eedge = edgePot(eidx);
if opt.DEBUG
    fprintf('Node = %.2f, Edge = %.2f, Sum = %.2f\n', sum(Enode), sum(Eedge), sum(Enode)+sum(Eedge));
    fprintf('Esmooth = %.2f\n', sum(Esmooths(eidx)));
end

%% update trace
newswc = swct;
newswc(:,5) = states(optimalDecoding);
end

