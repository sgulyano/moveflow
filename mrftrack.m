function [ swc, flag ] = mrftrack( swc, Iori, soma, soma_prev, swc_model, soma_model, opt )
% MRFTRACK
% By Sarun Gulyanon
%
% Usage
%
% [final_maskA, final_maskB] = main_coseg( IA, IB, WA, WB, maskA, maskB,...
%                               somaA, somaB, maskA_gold, maskB_gold, opt )
%
% This function tracking neuron in calcium image datasets using Markov 
% Network where dendrites are modeled as the articulated body. MRF is 
% optimized using fusion moves (QPBO + alpha expansion).
%
% Inputs:
%   swc           trace at time t_{n-1} in SWC format
%   Iori          calcium image stack at time t_n
%   soma          soma location at time t_n
%   soma_prev     soma location at time t_{n-1}
%   swc_model     trace at time t_0 in SWC format
%   soma_model    soma location at time t_0
%   opt      miscellaneous parameters.
%       'alpha_d'   Tune dynamic model. Default 1/5.
%       'alpha_t'   Tune transition similarity. Default 1/10.
%       'DEBUG'     Enable debugging mode. Default false.
%       'MAXDIST'   Maximum displacement. Default 10.
%
% Outputs:
%   swc           Trace at time t_n in SWC format
%   flag          True if the neuron moves out-of-bound, otherwise False.

if nargin < 7;  opt = struct();  end;
if ~isfield(opt,'alpha_d');       opt.alpha_d 	= 1/5;      end;
if ~isfield(opt,'alpha_t');       opt.alpha_t 	= 1/10;     end;

if ~isfield(opt,'DEBUG');       opt.DEBUG   = false;    end;
if ~isfield(opt,'MAXDIST');     opt.MAXDIST = 10;       end;

if all(soma == 0)
    flag = true;
    swc = [];
    return;
end

I = Iori;

%% init : find transition model
prev_pnt = swc(:,3:4);
numVar = size(prev_pnt,1);

% by soma location of given trace
dsoma_model = soma - soma_model;
model_pnt = swc_model(:,3:4) + ones(numVar,1)*dsoma_model(1:2);

% by soma location
dsoma = soma - soma_prev;
init_pnt = prev_pnt + ones(numVar,1)*dsoma(1:2);
pnt = round(init_pnt);

%% registration using MRF over tree structure
% Generate Adjacency Matrix
MAXDIST = opt.MAXDIST;
[xx, yy] = meshgrid(-MAXDIST:MAXDIST,-MAXDIST:MAXDIST);
states = [xx(:), yy(:)];
numLb = size(states,1);

% neuron graph
nb = swc(:,[1 7]);
nb = nb(all(nb>0,2),:);
deg = histc(nb(:),1:numVar);

mean_len = mean(sqrt(sum((swc(nb(:,1),3:4)-swc(nb(:,2),3:4)).^2,2)));
mean_len = max(mean_len, 3);

T = sparse(nb(:,1), nb(:,2), 1, numVar, numVar);
T = T + T';

% add shape smoothness edge to nodes' second level neighbor
nb1 = arrayfun(@(x)nchoosek(find(T(x,:)),2), find(deg>1), 'UniformOutput', false);
nb1len = cellfun(@(x)(size(x,1)*ones(size(x,1),1)), nb1, 'UniformOutput', false);
nb1 = vertcat(nb1{:});
nb1len = vertcat(nb1len{:});

G1 = sparse(nb1(:,1), nb1(:,2), nb1len, numVar, numVar); G1 = G1 + G1';

% find neuron branch id
root = find(deg ~= 2, 1);
[~,~,pred] = bfs(T, root);
brid = zeros(numVar,1);
num = 0;
for i = find(deg==1)'
    if i == root, continue, end;
    branch = path_from_pred(pred, i);
    brind = [0; find(deg(branch)>2); length(branch)];
    for j = 2:length(brind)
        if brid(branch(brind(j))) > 0, continue; end
        lst = branch(brind(j-1)+1:brind(j));
        if isempty(lst), continue, end;
        num = num+1;
        brid(lst) = num;
    end
end

% add repulsive edge at leaf nodes to nodes from diff branch
nb2 = [];
dist = zeros(numVar,numVar);
for i = 1:numVar
    dist(i,:) = sqrt(sum(bsxfun(@minus, pnt(i,:), pnt).^2,2));
end

% figure;
% xx = [swc(nb(:,1),3), swc(nb(:,2),3), nan(size(nb,1),1)]';
% yy = [swc(nb(:,1),4), swc(nb(:,2),4), nan(size(nb,1),1)]';
% plot(xx(:), yy(:), '-o', 'LineWidth', 2);
% set(gca,'Ydir','reverse')
% hold on;
for i = 1:num
    pi = find(brid == i & deg < 3);
    if isempty(pi), continue, end;
    for j = i+1:num
        pj = find(brid == j & deg < 3);
        if isempty(pj), continue, end;
        
        if length(pi) < length(pj)
            [dmini, posi] = min(dist(pi,pj),[],2);
            [dmin_sort, order] = sort(dmini);
            [~,ia,~] = unique(posi(order), 'first');
            
            ia( dmin_sort(order(ia)) > 2*MAXDIST ) = [];
            idx = order(ia);
            nb_ij = [pi(idx) pj(posi(idx))];
        else
            [dminj, posj] = min(dist(pi,pj),[],1);
            [dmin_sort, order] = sort(dminj);
            [~,ia,~] = unique(posj(order), 'first');
            
            ia( dmin_sort(order(ia)) > 2*MAXDIST ) = [];
            idx = order(ia);
            nb_ij = [pj(idx) pi(posj(idx))];
        end
        
        if isempty(nb_ij), continue; end
        nb_ij_ind = sub2ind(size(G1), nb_ij(:,1), nb_ij(:,2));
        nb_ij(full(G1(nb_ij_ind)>0),:) = [];
        
        if ~isempty(nb_ij)
%             xx = [swc(nb_ij(:,1),3), swc(nb_ij(:,2),3), nan(size(nb_ij,1),1)]';
%             yy = [swc(nb_ij(:,1),4), swc(nb_ij(:,2),4), nan(size(nb_ij,1),1)]';
%             plot(xx(:), yy(:),'LineWidth',1);
%             keyboard;
            nb2 = [nb2; nb_ij];
        end
    end
end
% hold off;

G2 = sparse(nb2(:,1), nb2(:,2), 1, numVar, numVar); G2 = max(G2, G2');

if full(any(any( T>0 & G2>0 )) | any(any( G1>0 & G2>0 )) | any(any( T>0 & G1>0 )))
    keyboard;
    error('edge overlap')
    
end
edgeStruct = UGM_makeEdgeStruct(T+G1+G2,numLb);

% create mask
swc2 = swc; swc2(:,3:4) = model_pnt;
mask = swc2pixel( swc2, [size(I) max(swc2(:,5))+1] );
mask = imdilate(mask, ones(2*MAXDIST+1));

Iadj = imadjust(I, [.02, .4], [0 1]);

% find vesselness
Options.BlackWhite = false;
Options.FrangiScaleRange = [1 3];
Options.FrangiScaleRatio = 0.5;
[Ivessel,~,~,Ibr]=FrangiFilter2D(Iadj*255,Options);
Iv = Ivessel;
Ivessel(~mask) = 0;
Ivessel = Ivessel ./ max(Ivessel(:));
Ibr(~mask) = 0;
Ibr = Ibr ./ max(Ibr(:));

% get linear index of nodes
p_ind = ones(numVar,numLb);
p_x = bsxfun(@plus, pnt(:,1), states(:,1)');
p_y = bsxfun(@plus, pnt(:,2), states(:,2)');
idxOut = p_x < 1 | p_x > size(I,2) | p_y < 1 | p_y > size(I,1);
p_ind(~idxOut) = sub2ind(size(I), p_y(~idxOut), p_x(~idxOut));

%% add unaries 
% Image feature: vesselness & function of eigenvalues of Hessian matrix
Eimg = zeros(numVar,numLb);
Eimg(deg>2,:) = Ibr(p_ind(deg>2,:));
Eimg(deg<=2,:) = Ivessel(p_ind(deg<=2,:));
Eimg(idxOut) = eps;

% Dynamics model (distance to model)
dx = bsxfun(@plus, pnt(:,1)-model_pnt(:,1), states(:,1)');
dy = bsxfun(@plus, pnt(:,2)-model_pnt(:,2), states(:,2)');
dist = sqrt(dx.^2+dy.^2);
Emodel = 1.5*exp(-dist * opt.alpha_d);

nodePot = Eimg + Emodel;

%% add pairwise potential
edgePot = ones(numLb,numLb,edgeStruct.nEdges);
if opt.DEBUG
    Eshapes = zeros(numLb,numLb,edgeStruct.nEdges);
    Esmooths = zeros(numLb,numLb,edgeStruct.nEdges);
    Erepul = zeros(numLb,numLb,edgeStruct.nEdges);
end
% Transition smoothness
vx = bsxfun(@minus, states(:,1), states(:,1)');
vy = bsxfun(@minus, states(:,2), states(:,2)');
dsmooth = sqrt(vx.^2+vy.^2);
Esmooth = exp(-dsmooth * opt.alpha_t); % ones(size(vx));
Esmooth = bsxfun( @rdivide, Esmooth, 2*mean(Esmooth(:)) );

for e = 1:edgeStruct.nEdges
    i = edgeStruct.edgeEnds(e,1);
    j = edgeStruct.edgeEnds(e,2);
    if T(i,j) == 1
        % Shape smoothness
        dx = bsxfun(@minus, pnt(i,1)+states(:,1), pnt(j,1)+states(:,1)');
        dy = bsxfun(@minus, pnt(i,2)+states(:,2), pnt(j,2)+states(:,2)');
        dshape = dx.^2+dy.^2;
        Eshape = -dshape/(MAXDIST^2);
        
        edgePot(:,:,e) = Esmooth + Eshape;
        if opt.DEBUG
            Eshapes(:,:,e) = Eshape;
            Esmooths(:,:,e) = Esmooth;
        end
    elseif G1(i,j) > 0
        % Shape smoothness
        dx = bsxfun(@minus, pnt(i,1)+states(:,1), pnt(j,1)+states(:,1)');
        dy = bsxfun(@minus, pnt(i,2)+states(:,2), pnt(j,2)+states(:,2)');
        dshape = dx.^2+dy.^2;
        Eshape = dshape/(G1(i,j)*4*MAXDIST^2);
        
        edgePot(:,:,e) = Eshape;
        if opt.DEBUG
            Eshapes(:,:,e) = Eshape;
        end
    elseif G2(i,j) == 1
        % Part-based model (repulsive edge)
        dx = bsxfun(@minus, pnt(i,1)+states(:,1), pnt(j,1)+states(:,1)');
        dy = bsxfun(@minus, pnt(i,2)+states(:,2), pnt(j,2)+states(:,2)');
        Erepulse = 10*double(sqrt(dx.^2+dy.^2) > mean_len); % ones(size(vx));
        
        edgePot(:,:,e) = Erepulse;
        if opt.DEBUG
            Erepul(:,:,e) = Erepulse;
        end
    else
        error('unknown edge type');
    end
end

%% inference
% optimalDecoding = UGM_Decode_Tree(nodePot,edgePot,edgeStruct);
% optimalDecoding = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
[~, optimalDecoding] = fusion_moves(-nodePot, -edgePot, edgeStruct.edgeEnds, 20);

nidx = sub2ind(size(nodePot),1:numVar,optimalDecoding');
Enode = nodePot(nidx);
eidx = sub2ind(size(edgePot),...
        optimalDecoding(edgeStruct.edgeEnds(:,1))',...
        optimalDecoding(edgeStruct.edgeEnds(:,2))',...
        1:edgeStruct.nEdges);
Eedge = edgePot(eidx);
fprintf('Node = %.2f, Edge = %.2f, Sum = %.2f\n', sum(Enode), sum(Eedge), sum(Enode)+sum(Eedge));
fprintf('Eimg = %.2f, Edist = %.2f\n', sum(Eimg(nidx)), sum(Emodel(nidx)));
if opt.DEBUG
    fprintf('Eshape = %.2f, Esmooth = %.2f, Erepul = %.2f\n', sum(Eshapes(eidx)), sum(Esmooths(eidx)), sum(Erepul(eidx)));
end
%% plot
if opt.DEBUG
    [ii, jj] = find(mask);
    xmin = max(min(jj)-MAXDIST, 1); ymin = max(min(ii)-MAXDIST, 1);
    xmax = min(max(jj)+MAXDIST, size(I,2)); ymax = min(max(ii)+MAXDIST, size(I,1));
    Ivcrop = Iv(ymin:ymax, xmin:xmax);
    Icrop = Iori(ymin:ymax, xmin:xmax);

    swc2 = swc;
    [~, nodeDecoding] = max(nodePot, [], 2);

    figure(7); subplot(2,1,1); imshow(Iadj); title('adjust intensity image');
    rectangle('Position',[xmin,ymin,xmax-xmin,ymax-ymin],...
             'LineWidth',2,'LineStyle','--','EdgeColor','g')
    viscircles(soma(1:2), soma(3), 'EdgeColor', 'r');

    swc2(:,3:4) = bsxfun(@minus, init_pnt + states(nodeDecoding,:), [xmin ymin]-1);
    subplot(2,2,3); EV_plot_img( Ivcrop, swc2, 'c-o' ); title('Only unary(R=model,Ives)');
    swc2(:,3:4) = bsxfun(@minus, model_pnt, [xmin ymin]-1);
    hold on; EV_plot_img( [], swc2, 'r--' );

    swc2(:,3:4) = bsxfun(@minus, init_pnt + states(optimalDecoding,:), [xmin ymin]-1);
    subplot(2,2,4); EV_plot_img( Icrop, swc2, 'c-o' ); title('MRF (R=init,Is)');
    swc2(:,3:4) = bsxfun(@minus, init_pnt, [xmin ymin]-1);
    hold on; EV_plot_img( [], swc2, 'r--' ); 
end

%% update trace
swc(:,3:4) = init_pnt + states(optimalDecoding,:);
% check if half of neuron is out-of-bound, then tracking is done
s = size(I);
maxsize = ones(numVar,1)*s([2 1]);
count_out = sum(any(pnt < 1 | pnt > maxsize, 2));
if count_out > numVar/2
    flag = true;
else
    flag = false;
end
end
