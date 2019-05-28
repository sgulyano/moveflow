function cursnk = AC_remesh(cursnk)
% AC_REMESH     Remesh an active contour (AC) to ensure the resolution.
%     vertex1 = AC_REMESH(vertex0)
%     vertex1 = AC_REMESH(vertex0,res)
%     vertex1 = AC_REMESH(vertex0,res,method)
%     vertex1 = AC_REMESH(vertex0,res,method,type)
% 
%     Inputs
%     cursnk      a struct of snake consisting of vert, kstr, iter
%     vertex0     position of the vertices, n-by-2 matrix, each row of 
%                 which is [x y]. n is the number of vertices.
%     res         desired resolution (distance between vertices) in pixels,
%                 default value is 1.
%     method      'adaptive' - (default) the distances between vertices
%                           will be larger than res/2 and less than res*2,
%                           the shape will be the same
%                 'equal' - the distances between vertices will be the same, 
%                           the shape may be slightly different after remeshing
%     type        'close' - close contour (default), the last vertex and
%                           first vertex are connected 
%                 'open'  - open contour,  the last vertex and first vertex
%                           are not connected 
%     Outputs
%     vertex1     position of the vertices after remeshing, m-by-2 matrix
%     
%     Example
%         See EXAMPLE_VFC, EXAMPLE_PIG.
%
%     See also AMT, AC_DEFORM, AM_VFC, AM_VFK, AM_PIG, AC_INITIAL, 
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
%   02-18-2006  support open type
%   01-30-2009  minor bug fix

%% inputs check
% if ~ismember(nargin, 3:5) || numel(res) ~= 1,
%     error('Invalid inputs to AC_REMESH!')    
% end
if size(cursnk.vert,2) ~= 3 
    error('Invalid vertex matrix!')
end

%% delete the vertices close to each other (distance less than res/2)
N = size(cursnk.vert,1);
if N < 3
    return
end

%% remesh
pnt = cursnk.vert;
dif = diff(pnt);
pnt(sum(dif,2)==0,:) = [];  % Ensure  strict monotonic order.
dif(sum(dif,2)==0,:) = [];  % consistent to cursnk.vert
d = sqrt(sum(dif.^2,2));  % distance between adjacent points
D = sum(d);
di = linspace(0,D,round(D));
dcum = [0;cumsum(d)];
vert = zeros(round(D),3);
for i=1:3,
    vert(:,i) = interp1(dcum,pnt(:,i),di)';
end
cursnk.vert = vert;

end