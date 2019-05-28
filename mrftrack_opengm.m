function [ swc, flag ] = mrftrack_opengm( swc, Iori, soma, soma_prev, swc_model, soma_model, opt  )
%MRFTRACK tracking neuron using 3rd order Markov Network with cycle

if nargin < 7;  opt = struct();  end;
if ~isfield(opt,'alpha');       opt.alpha     = 1.0;    end;
if ~isfield(opt,'beta');        opt.beta      = 1.0;    end;
if ~isfield(opt,'oof_thr');     opt.oof_thr   = 2;      end;
if ~isfield(opt,'r_oof');       opt.r_oof     = 4;      end;
if ~isfield(opt,'sigmas');      opt.sigmas    = 2;      end;

if ~isfield(opt,'OF_thr');      opt.OF_thr = 50;        end;
if ~isfield(opt,'DEBUG');       opt.DEBUG     = false;  end;

if ~isfield(opt,'python'); opt.python = '/Users/Yo/anaconda2/bin/python'; end;

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
MAXDIST = 10;
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

% find neuron branch id
[~,~,pred] = bfs(T, 1);
brid = zeros(numVar,1);
num = 0;
for i = find(deg==1)'
    if i == 1, continue, end;
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
for i = find(deg==1)'
    ds = sqrt(sum(bsxfun(@minus, pnt(i,:), pnt).^2,2));
    pnb = full(T(:,i)); pnb(i) = 1;
    cand = find(ds < MAXDIST & ~pnb & brid~=brid(i));
    nb2 = [nb2; i*ones(size(cand)) cand];
end

G2 = sparse(nb2(:,1), nb2(:,2), 1, numVar, numVar); G2 = max(G2, G2');
if full(any(any( T>0 & G2>0 )))
    error('edge overlap')
end
edgeStruct = UGM_makeEdgeStruct(T+G2,numLb);

% create mask
swc2 = swc; swc2(:,3:4) = model_pnt;
mask = swc2pixel( swc2, [size(I) max(swc2(:,5))+1] );
mask = imdilate(mask, ones(2*MAXDIST+1));

Iadj = imadjust(I, [.02, .4], [0 1]);
figure(11), imshow(Iadj)

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


%% init graphical model
[~, executable, isloaded] = pyversion;
if ~strcmp(executable, opt.python) && ~isloaded
    pyversion(opt.python)
end

gm = py.opengm.gm( numLb*ones(1,numVar,'int32') );

%% add unaries (pixel intensity of smooth img = intensity around pnt)
Eimg = zeros(numVar,numLb);
Eimg(deg>2,:) = Ibr(p_ind(deg>2,:));
Eimg(deg<=2,:) = Ivessel(p_ind(deg<=2,:));
Eimg(idxOut) = eps;

% local smooth transition (distance to model)
dx = bsxfun(@plus, pnt(:,1)-model_pnt(:,1), states(:,1)');
dy = bsxfun(@plus, pnt(:,2)-model_pnt(:,2), states(:,2)');
dist = sqrt(dx.^2+dy.^2);
Edist = exp(-dist/5);

% node potential
nodePot = -(Eimg + Edist);

% add unaries
nodePotT = nodePot';
unaryFunctionIds = gm.addFunctions( py.numpy.array(nodePotT(:)').reshape(numVar,numLb) );
gm.addFactors(unaryFunctionIds,py.numpy.arange(numVar));

%% add edge potential
% global smooth transition
vx = bsxfun(@minus, states(:,1), states(:,1)');
vy = bsxfun(@minus, states(:,2), states(:,2)');
dsmooth = sqrt(vx.^2+vy.^2);
% Esmooth = ones(size(vx));
Esmooth = exp(-dsmooth/10);
Esmooth = -bsxfun(@rdivide, Esmooth, 2*mean(Esmooth(:)) );
% EsmoothT = Esmooth';
% EsmoothId = gm.addFunction( py.numpy.array(EsmoothT(:)').reshape(numLb,numLb) );

for e = 1:edgeStruct.nEdges
    i = edgeStruct.edgeEnds(e,1);
    j = edgeStruct.edgeEnds(e,2);
    vis = [i-1 j-1];
    if T(i,j) == 1
        dx = bsxfun(@minus, pnt(i,1)+states(:,1), pnt(j,1)+states(:,1)');
        dy = bsxfun(@minus, pnt(i,2)+states(:,2), pnt(j,2)+states(:,2)');
        dshape = sqrt(dx.^2+dy.^2);
        Eshape = exp(-dshape/10);
        
        edgePot = Esmooth + Eshape;
        edgePot = -bsxfun(@rdivide, edgePot, 2*mean(edgePot(:)) );
        edgePotT = edgePot';
        
        edgePotId = gm.addFunction( py.numpy.array(edgePotT(:)').reshape(numLb,numLb) );
        gm.addFactor(edgePotId,int32(vis));
        
%         gm.addFactor(EsmoothId,int32(vis));
    elseif G2(i,j) == 1
        % repulsive edge
        dx = bsxfun(@minus, pnt(i,1)+states(:,1), pnt(j,1)+states(:,1)');
        dy = bsxfun(@minus, pnt(i,2)+states(:,2), pnt(j,2)+states(:,2)');
%         Erepulse = ones(size(vx));
        Erepulse = -double(sqrt(dx.^2+dy.^2) > mean_len);
        ErepulseT = Erepulse';
        
        ErepulseId = gm.addFunction( py.numpy.array(ErepulseT(:)').reshape(numLb,numLb) );
        gm.addFactor(ErepulseId,int32(vis));
    else
        error('unknown edge type');
    end
end

% %% add second order potential (local shape smoothness)
% tic;
% for q = find(deg>1)'
%     nbs = nchoosek(find(T(q,:)),2);
%     for j = 1:size(nbs,1)
%         p = nbs(j,1);
%         r = nbs(j,2);
%         dpqx = bsxfun(@minus, pnt(p,1)+states(:,1), pnt(q,1)+states(:,1)');
%         dpqy = bsxfun(@minus, pnt(p,2)+states(:,2), pnt(q,2)+states(:,2)');
%         dqrx = bsxfun(@minus, pnt(q,1)+states(:,1), pnt(r,1)+states(:,1)');
%         dqry = bsxfun(@minus, pnt(q,2)+states(:,2), pnt(r,2)+states(:,2)');
%         
%         d_x = repmat(dpqx, [1 1 numLb]) - repmat(reshape(dqrx, [1 numLb numLb]), [numLb 1 1]);
%         d_y = repmat(dpqy, [1 1 numLb]) - repmat(reshape(dqry, [1 numLb numLb]), [numLb 1 1]);
%         
%         Eshape = -double(sqrt(d_x.^2+d_y.^2)) / 5;
%         EshapeT = permute(Eshape, [2 1 3]);
%         
%         
%         vis = [nbs(j,1) i nbs(j,2)];
%         [vis, order] = sort(vis);
%         EshapeT = permute(EshapeT, order);
%         EshapeId = gm.addFunction( py.numpy.array(EshapeT(:)').reshape(numLb,numLb,numLb) );
%         gm.addFactor(EshapeId,int32(vis));
%     end
% end
% toc;

% nb1 = arrayfun(@(x)(nchoosek(find(T(x,:)),2)), nonleaf, 'UniformOutput', false);
% % shape_edge_num = cellfun(@(x)size(x,1), nb1);
% nb1 = vertcat(nb1{:});
% 
% % shape_edge_n = sparse(nb1(:,1), nb1(:,2), shape_edge_num, numVar, numVar);


%% inference using BP
gm
parameter = py.opengm.InfParam(pyargs('steps',100,'damping',0.9,'convergenceBound',0.001));
inf2 = py.opengm.inference.TreeReweightedBp(gm,pyargs('parameter',parameter));
% parameter = py.opengm.InfParam(pyargs('rounds',10,'useKovtunsMethod',false,'strongPersistency',true));
% inf2 = py.opengm.inference.Mqpbo(gm,pyargs('parameter',parameter));
inf2.infer();
decoding = inf2.arg();

optimalDecoding = double( py.array.array('d', py.numpy.nditer(decoding)) );  % Add order='F' to get data in column-major order (as in Fortran 'F' and Matlab

%%
nidx = sub2ind(size(nodePot),1:numVar,optimalDecoding+1);
Enode = nodePot(nidx);
fprintf('Node = %f : Eimg = %d, Edist = %f\n', sum(Enode), sum(Eimg(nidx)), sum(Edist(nidx)));

%%
[ii, jj] = find(mask);
xmin = max(min(jj)-MAXDIST, 1); ymin = max(min(ii)-MAXDIST, 1);
xmax = min(max(jj)+MAXDIST, size(I,2)); ymax = min(max(ii)+MAXDIST, size(I,1));
Ivcrop = Iv(ymin:ymax, xmin:xmax);
Icrop = Iori(ymin:ymax, xmin:xmax);

swc2 = swc;
[~, nodeDecoding] = min(nodePot, [], 2);
swc2(:,3:4) = bsxfun(@minus, init_pnt + states(nodeDecoding,:), [xmin ymin]-1);
figure(2); subplot(1,2,1); EV_plot_img( Ivcrop, swc2, 'c-o' ); title('Unary potential only');
% swc2(:,3:4) = bsxfun(@minus, init_pnt, [xmin ymin]-1);
% hold on; EV_plot_img( [], swc2, 'r--o' );
swc2(:,3:4) = bsxfun(@minus, model_pnt, [xmin ymin]-1);
hold on; EV_plot_img( [], swc2, 'r--o' );

swc2(:,3:4) = bsxfun(@minus, init_pnt + states(optimalDecoding+1,:), [xmin ymin]-1);
subplot(1,2,2); EV_plot_img( Icrop, swc2, 'c-o' ); title('MRF');
swc2(:,3:4) = bsxfun(@minus, init_pnt, [xmin ymin]-1);
hold on; EV_plot_img( [], swc2, 'r--' ); 
% swc2(:,3:4) = bsxfun(@minus, model_pnt, [xmin ymin]-1);
% hold on; EV_plot_img( [], swc2, 'g--' ); 
%%
% % hold on; EV_plot_img( [], swc_mod, 'g--' ); 

swc(:,3:4) = init_pnt + states(optimalDecoding+1,:);
% swc(:,3:4) = init_pnt;
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
