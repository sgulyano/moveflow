function cursnk = AC_deform(cursnk,alpha,beta,tau,Fext,...
        ev,mask,tube,snkcode,model,type)
% AC_DEFORM     Deform an active contour (AC), also known as snake.
%     vertex1 = AC_DEFORM(vertex0,alpha,beta,tau,Fext,iter)
%     vertex1 = AC_DEFORM(vertex0,alpha,beta,tau,Fext,iter,type)
%   
%     Inputs
%     cursnk      a struct of a snake
%     alpha       AC elasticity (1st order) parameter ranges from 0 to 1.
%     beta        AC rigidity (2nd order) parameter ranges from 0 to 1.
%     tau         time step of each iteration.
%     Fext        the external force field,d1-by-d2-by-2 matrix, 
%                 the force at (x,y) is [Fext(y,x,1) Fext(y,x,2)].
%     ev          eigenvector at each pixels; d1-by-d2-by-d3-by-3 matrix,
%     mask        binary image of neuron
%     tube        probability map of being neuron
%     snkcode     logical value indicate snake's end collision
%     model       force from Neuron model
%     type        'close' - close contour (default), the last vertex and
%                           first vertex are connected 
%                 'open'  - open contour,  the last vertex and first vertex
%                           are not connected 
%               
%     Outputs
%     cursnk      snake after deformation, n-by-2 matrix
%     
%     Note that if the vertices are outside the valid range, i.e., y>d1 ||
%     y<1 || x>d2 || x<1, they will be pulled inside the valid range. 
% 
%     Example
%         See EXAMPLE_VFC, EXAMPLE_PIG.
%
%     See also AMT, AM_VFC, AM_VFK, AM_PIG, AC_INITIAL, AC_REMESH,
%     AC_DISPLAY, AM_GVF, EXAMPLE_VFC, EXAMPLE_PIG. 
% 
%     Reference
%     [1] Bing Li and Scott T. Acton, "Active contour external force using
%     vector field convolution for image segmentation," Image Processing,
%     IEEE Trans. on, vol. 16, pp. 2096-2106, 2007.  
%     [2] Bing Li and Scott T. Acton, "Automatic Active Model
%     Initialization via Poisson Inverse Gradient," Image Processing,
%     IEEE Trans. on, vol. 17, pp. 1406-1420, 2008.   
% 
% (c) Copyright Bing Li 2005 - 2009.

% Revision Log
%   11-30-2005  original 
%   01-30-2006  external force interpolation outside the image 
%   02-18-2006  add open contour codes
%   01-30-2009  minor bug fix

% inputs check
if nargin == 6
    type = 'close';
end

N = size(cursnk.vert,1);
if size(cursnk.vert,2) ~= 3
    error('Invalid vertex matrix!')
end

if N < 3
    return
end

%% compute T = (I + tao*A) of equation (9) in reference [1]
Lap = sparse(1:N, 1:N, -2) + sparse(1:N, [N 1:N-1], 1) + sparse(1:N, [2:N 1], 1);

if strcmp(type,'open'), % offset tau for boundary vertices    
    tau = sparse(1:N, 1:N, tau);
    tau(1) = 0;
    tau(end) = 0;
    offset = sparse(1:N, 1:N, 1);
    offset(1,1)=0;      offset(N,N) = 0;
    offset(1,2)=1;      offset(N,N-1) = 1;
    Lap = offset*Lap;
end

T =sparse(1:N,1:N,1)+ tau*(beta*Lap*Lap-alpha*Lap);
% hack
tau(1) = 0.5;
tau(end) = 0.5;

%% Another way to compute T for close AC
% a = beta;
% b = -alpha - 4*beta;
% c = 2*alpha + 6*beta;
% 
% T = sparse(1:N,1:N,1) + tau*(sparse(1:N,1:N,c) + ...
%     sparse(1:N,[N,1:N-1],b) + sparse(1:N,[2:N,1],b)...
%     + sparse(1:N,[N-1,N,1:N-2],a) + sparse(1:N,[3:N,1,2],a));

%% Deform
center = size(tube)/2;
center = center([2,1,3]);
F = zeros(size(cursnk.vert)); % image energy (GVF) term

% check if snake points lie outside the image
IdxOut = cursnk.vert(:,1)<1 | cursnk.vert(:,1)>size(tube,2) ...
        | cursnk.vert(:,2)<1 | cursnk.vert(:,2)>size(tube,1) ...
        | cursnk.vert(:,3)<1 | cursnk.vert(:,3)>size(tube,3);
IdxIn = ~IdxOut;
F(IdxOut,:) = 0;
% get image energy (GVF) term
for i=1:3	
	% interpolate the external force for vertices within the range
	F(IdxIn,i) = ba_interp3(Fext{i},...
            cursnk.vert(IdxIn,1),cursnk.vert(IdxIn,2),cursnk.vert(IdxIn,3)...
            ,'linear');
	% for points outside image, pointing back to volume
	Idxi = cursnk.vert(:,i)<1;
	F(Idxi,i) = 0;
	Idxi = cursnk.vert(:,i)>(center(i)*2);
	F(Idxi,i) = 0;
end
if ~isempty(IdxOut)
	% normalize the forces outside the image
    F(IdxOut,:) = 0;
end       
% get stretching direction
cs1 = -(cursnk.vert(2,:) - cursnk.vert(1,:));
csN = cursnk.vert(N,:) - cursnk.vert(N-1,:);

ev1N = zeros(2,3);
for i=1:3
    ev1N(:,i) = ba_interp3(ev{i}, ...
            [cursnk.vert(1,1), cursnk.vert(N,1)],...
            [cursnk.vert(1,2), cursnk.vert(N,2)],...
            [cursnk.vert(1,3), cursnk.vert(N,3)],'linear');
end
Estr = zeros(size(F));
Estr(1,:) = sign(cs1*ev1N(1,:)')*ev1N(1,:)/norm(ev1N(1,:));
Estr(N,:) = sign(csN*ev1N(2,:)')*ev1N(2,:)/norm(ev1N(2,:));
Estr(isnan(Estr)) = 0;


% get stretching coefficient
kstr = zeros(size(Estr));
p1 = cursnk.vert(1,:);
p1r = round(p1);
% using mask for leakage prevention and expS for stretch confidence
try
    kstr(1,:) = 0.5 * double(mask(p1r(2),p1r(1),p1r(3))) * ...
            (1 + ba_interp3(tube,p1(1),p1(2),p1(3),'linear'));
catch
    kstr(1,:) = 0;
end

pN = round(cursnk.vert(N,:));
pNr = round(pN);
try
    kstr(N,:) = 0.5 * double(mask(pNr(2),pNr(1),pNr(3))) * ...
            (1 + ba_interp3(tube,pN(1),pN(2),pN(3),'linear'));
catch
    kstr(N,:) = 0;
end
kstr(isnan(kstr)) = 0;

% get external force coefficient
k = ones(size(Estr));
if ~snkcode(1) || cursnk.length > 1.5*N
    k(1,:) = 0;
end

if ~snkcode(2) || cursnk.length > 1.5*N
    k(N,:) = 0;
end

% get model force
Emod = zeros(size(F));
Emod(1,:) = model.confidenceH .* model.predict_directionH;
Emod(N,:) = model.confidenceT .* model.predict_directionT;
% Emod(1,:) = 0.5 * double(mask(p1r(2),p1r(1),p1r(3))) * model.confidenceH .* model.predict_directionH;
% Emod(N,:) = 0.5 * double(mask(pNr(2),pNr(1),pNr(3))) * model.confidenceT .* model.predict_directionT;

% External force
% equation (9) in reference [1]
dEext = k.*(F+(kstr.*Estr)+Emod);

% update snake vertices
cursnk.vert = T\(cursnk.vert+tau*dEext);
