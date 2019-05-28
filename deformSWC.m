function [ swcs_ddad, swcs_ddae ] = deformSWC( V, swcs, track_time, soma_pos, configs, spline, opt  )
%DEFORMSWC deform SWC by our neuron's plane model in 2D
if nargin < 7;  opt = struct();  end;
if ~isfield(opt,'dx');          opt.dx = 4;                         end;

swcs_ddad = cell(1,size(V,3));
swcs_ddae = cell(1,size(V,3));

dx = opt.dx;
sizeI = [size(V,1), size(V,2)];

% locate soma
ddad_somaY = round(swcs{1}(1,4));
ddae_somaY = round(swcs{2}(1,4));
somaY = round((ddad_somaY + ddae_somaY)/2);

% create baseline mask
I_mask1 = swc2pixel( swcs{1}, sizeI );
I_mask2 = swc2pixel( swcs{2}, sizeI );
I2D_mask = imdilate(I_mask1 | I_mask2, ones(3));

% get bounding box
[ii, jj] = find(I2D_mask);
xmin = min(jj) - 1;
yrange = [min(ii) max(ii)];
ymin = yrange(1) - 1;

swc_ddad = swcs{1}; swc_ddad(:,3:4) = round(swc_ddad(:,3:4));
swcs_ddad{track_time+1} = swc_ddad;

swc_ddae = swcs{2}; swc_ddae(:,3:4) = round(swc_ddae(:,3:4));
swcs_ddae{track_time+1} = swc_ddae;


%% create animation
for t = track_time+2:size(V,3)
    if isempty(configs{t,1})
        break;
    end
    ddadX = soma_pos{1}(t,1);
    ddadY = soma_pos{1}(t,2);
    ddaeX = soma_pos{2}(t,1);
    ddaeY = soma_pos{2}(t,2);
    newsomaX = (ddadX + ddaeX)/2;
    dsomaY = round(somaY - (ddadY+ddaeY)/2);
    
    Ox = ffd_init_from_config(configs(t,:), newsomaX, dx);
    
    % transform dendrite in 2D
    Tx = ffd_interpolate(Ox, spline);
    newyrange = yrange - dsomaY;
    [cx, cy] = meshgrid(round(Tx), newyrange(1):newyrange(2));
    
    swc_ddad = swcs{1};
    swc_ddad(:,3) = round( interp2(cx, swcs{1}(:,3) - xmin, swcs{1}(:,4) - ymin) );
    swc_ddad(:,4) = round( interp2(cy, swcs{1}(:,3) - xmin, swcs{1}(:,4) - ymin) );
    swcs_ddad{t} = swc_ddad;
    
    swc_ddae = swcs{2};
    swc_ddae(:,3) = round( interp2(cx, swcs{2}(:,3) - xmin, swcs{2}(:,4) - ymin) );
    swc_ddae(:,4) = round( interp2(cy, swcs{2}(:,3) - xmin, swcs{2}(:,4) - ymin) );
    swcs_ddae{t} = swc_ddae;
    
    figure(1); EV_plot_img(V(:,:,t), swcs_ddad{t});
    hold on; EV_plot_img([], swcs_ddae{t}, 'r'); hold off;
    title(t-1);
    pause(0.5);
end


end

