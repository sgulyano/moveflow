function snk = snake_deform(snk,alpha,beta,tau,Fext,ev,mask)
% SNAKE_DEFORM Deform an active contour (AC), also known as snake.
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
N = size(snk.vert,1);
if size(snk.vert,2) ~= 3
    error('Invalid vertex matrix!')
end

if N < 3
    return
end

%% compute T = (I + tao*A) of equation (9) in reference [1]
Lap = sparse(1:N, 1:N, -2) + sparse(1:N, [N 1:N-1], 1) + sparse(1:N, [2:N 1], 1);

% offset tau for boundary vertices
tau_tmp = tau;
tau = sparse(1:N, 1:N, tau);
tau(1) = 0;
tau(end) = 0;
offset = sparse(1:N, 1:N, 1);
offset(1,1)=0;      offset(N,N) = 0;
offset(1,2)=1;      offset(N,N-1) = 1;
Lap = offset*Lap;

T =sparse(1:N,1:N,1)+ tau*(beta*Lap*Lap-alpha*Lap);
% hack
tau(1) = tau_tmp;
tau(end) = tau_tmp; %tau_tmp; % * double(interp2(mask,snk.vert(end,1),snk.vert(end,2),'linear')<.5);

%% Another way to compute T for close AC
% a = beta;
% b = -alpha - 4*beta;
% c = 2*alpha + 6*beta;
% 
% T = sparse(1:N,1:N,1) + tau*(sparse(1:N,1:N,c) + ...
%     sparse(1:N,[N,1:N-1],b) + sparse(1:N,[2:N,1],b)...
%     + sparse(1:N,[N-1,N,1:N-2],a) + sparse(1:N,[3:N,1,2],a));

%% Deform
F = zeros(size(snk.vert)); % image energy (GVF) term

% check if snake points lie inside the image and mask
s = [size(Fext{1}) inf];
maxsize = ones(N,1)*s([2 1 3]);
IdxIn = (~any(snk.vert<1 | snk.vert>maxsize,2));% & interp2(double(mask),snk.vert(:,1),snk.vert(:,2),'linear',0)>.5;
% IdxIn([1 end]) = false;
% get image energy (GVF) term
for i=1:2
	% interpolate the external force for vertices within the range
	F(IdxIn,i) = interp2(Fext{i},snk.vert(IdxIn,1),snk.vert(IdxIn,2),'linear',0);
end

% % get stretching direction
% cs1 = -(snk.vert(2,:) - snk.vert(1,:));
% csN = snk.vert(N,:) - snk.vert(N-1,:);
% 
% ev1N = zeros(2,3);
% for i=1:2
%     ev1N(:,i) = interp2(ev{i},snk.vert([1 N],1),snk.vert([1 N],2),'linear');
% end
% Estr = zeros(size(F));
% Estr(1,:) = sign(cs1*ev1N(1,:)')*ev1N(1,:)/norm(ev1N(1,:));
% Estr(N,:) = sign(csN*ev1N(2,:)')*ev1N(2,:)/norm(ev1N(2,:));
% Estr(isnan(Estr)) = 0;
% 
% % get stretching coefficient
% kstr = zeros(size(Estr));
% pN = round(snk.vert(N,:));
% pNr = round(pN);
% % kstr(N,:) = double( mask(pNr(2),pNr(1),pNr(3)) );
% kstr(isnan(kstr)) = 0;

% External force
% equation (9) in reference [1]
dEext = F;%(F+(kstr.*Estr));

% update snake vertices
snk.vert = T\(snk.vert+tau*dEext);
% snk.vert(IdxIn,1:2) = vert(IdxIn,1:2);

% snk.vert(snk.vert < 1) = 1;
% snk.vert(snk.vert > maxsize) = maxsize(snk.vert > maxsize);
end