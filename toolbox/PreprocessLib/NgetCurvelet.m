
function [Ix,Iy,orientation] = NgetCurvelet(img,sig1)%,C,E)

% Don't forget that you can choose between diagonal real, diagonal complex
% and block complex thersholding by uncommenting portions of the code

% filename = 'C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\test.tif';

% filename = 'C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\hippocampal_006.tif';


% close all


%     out_filename = [filename(1:end-4)  '_out.tif'];
%    out_filename_cos = [filename(1:end-4) '_cos_out.tif'];
%    out_filename_sin = [filename(1:end-4) '_sin_out.tif'];
%    delete(out_filename);
%    delete(out_filename_cos);
%    delete(out_filename_sin);
    
%im_temp = squeeze(input(:,:,outer_counter));
%noisy_img = zeros(2^ceil(log2(size(im_temp,1))),2^ceil(log2(size(im_temp,2))));
%noisy_img(1:size(im_temp,1),1:size(im_temp,2)) = im_temp;
% noisy_img = double(imread(filename,40));

[M, N] = size(img);
noisy_img=img;
if(~exist('sig1','var'))
    sigma = 0.4*255;            % noise stdev is 10%
else
    sigma = sig1;
end

% Tuning parameters
finest = 1;                 % 1: curvelets at finest scale, 2: wavelets at finest scale
nbscales = log2(min(N,M)) - 3;
nbangles_coarse = 8;
nsigmas_coarse = 2.2;       % threshold proportional to nsigmas_coarse*sigma at all scales except finest
nsigmas_fine = 2.5;         % threshold proportional to nsigmas_fine*sigma at finest scale
nshifts = 1;                % number of translations (per dimension) considered in the cycle-spinning
nell1 = 0;                  % number of ell-1 iterations
neighb_weight = 0.5;        % for group thresholding, weight assigned to neighboring curvelets
tuning_neighb = 0.6;        % for group thresholding, parameter used to renormalize the weighted sum of coefficients squared
is_real = 0;

% get L2norm, put into E 
%F = ones(N,N);
F = ones(M,N); 
X = fftshift(ifft2(F)) * sqrt(numel(F));
disp('Computing L^2 norms ...');
tic;
C = fdct_wrapping(X,is_real,finest,nbscales,nbangles_coarse);

E = cell(size(C));
E2 = cell(size(C));
for s=1:length(C)
  E{s} = cell(size(C{s}));
  E2{s} = cell(size(C{s}));
  for w=1:length(C{s})
    A = C{s}{w};
    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / numel(A));
    E2{s}{w} = sqrt(A.*conj(A));
  end
end

disp('Computing parameters ...');
%[X_rows, X_cols, F_rows, F_cols, N_rows, N_cols] = fdct_wrapping_param(C,N,N);
[X_rows, X_cols, F_rows, F_cols, N_rows, N_cols] = fdct_wrapping_param(C,M,N);

% Cycle spinning
n = 0;
restored_img = 0;
for xshift = 1:nshifts
  for yshift = 1:nshifts
    
    shift_img = circshift(noisy_img,[xshift yshift]);
    n = n + 1;
    
    disp(['Direct transform, shift nr. ',num2str(n),' ...']);
    C = NfdctWrapping(shift_img,is_real,finest,nbscales,nbangles_coarse);
    Ccos = C;
    Csin = C;
    % Thresholding
    disp('Thresholding ...')
    thresh = nsigmas_coarse * sigma;
    for j = 1:length(C)
      if j == length(C), thresh = nsigmas_fine * sigma; end;
      for l = 1:length(C{j})
        
        thresh_jl = thresh*E{j}{l};
        
        % Uncomment for 'diagonal real thresholding'
        %                recjl = real(C{j}{l});     
        %                imcjl = imag(C{j}{l});
        %                recjl = recjl .* (abs(recjl) > thresh_jl);
        %                imcjl = imcjl .* (abs(imcjl) > thresh_jl);
        %                C{j}{l} = recjl + sqrt(-1)*imcjl;                
        
        % Uncomment for 'diagonal complex thresholding'
        %                modcjl = abs(C{j}{l});     
        %                argcjl = C{j}{l} ./ modcjl;
        %                modcjl = modcjl .* (modcjl > thresh_jl);
        %                C{j}{l} = argcjl .* modcjl;
        
        % Uncomment for 'block complex thresholding'
        modcjl = abs(C{j}{l});          
        argcjl = C{j}{l} ./ modcjl;
        rowstep = M/N_rows{j}{l};
        %colstep = N/N_cols{j}{l};
        colstep = N/N_cols{j}{l};
        evenquad = ~mod(ceil(l*4/length(C{j})),2);
        if evenquad,
          if (j == 1)||(finest==2 && j==nbscales), fcolsjl = 1; else fcolsjl = F_cols{j}{l}; end;
          rowshift = - round(F_rows{j}{l}/fcolsjl * rowstep);
          %[rowshift j l]
          %pause;
          testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[1 0]).^2 + circshift(modcjl,[-1 0]).^2 + ...
                                                    circshift(modcjl,[rowshift 1]).^2 + circshift(modcjl, [-rowshift -1]).^2));
        else
          if (j == 1)||(finest==2 && j==nbscales), frowsjl = 1; else frowsjl = F_rows{j}{l}; end;
          colshift =  round(-(F_cols{j}{l}/frowsjl * colstep));
          colshift=cast(colshift,'double');
          [colshift j l];
          circshiftVec=[1,colshift];
          %pause;
          debugcjl1=circshift(modcjl,[0 1]).^2;
          debugcjl2=circshift(modcjl,[0 -1]).^2;
          debugcjl3=circshift(modcjl,circshiftVec).^2;
          debugcjl4=circshift(modcjl,-circshiftVec).^2;
          testcjl=sqrt(modcjl.^2 + neighb_weight*(debugcjl1+debugcjl2+debugcjl3+debugcjl4));
          %testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[0 1]).^2 + circshift(modcjl,[0 -1]).^2 + ...
%                                                     circshift(modcjl,[1 -35]).^2 + circshift(modcjl, [-1 35]).^2));
        end;
        testcjl = testcjl ./ sqrt(1+4*neighb_weight*tuning_neighb);
%         figure,
%         subplot(2,1,1);imagesc(modcjl);
        modcjl = modcjl .* (testcjl > thresh_jl);
%         subplot(2,1,2);imagesc(modcjl);
        C{j}{l} = argcjl .* modcjl;
%         if j == 2
%            if l == outer_angles
%             C{j}{l} = C{j}{l};
%            else
%             C{j}{l} = 0*C{j}{l};
%            end
%         else
%             C{j}{l} = 0*C{j}{l};
%         end
%           if j ~= 1
%               C{j}{l} = 0*C{j}{l};
%           end
        if j ~=1 
            theta = (pi/4 -2*pi/numel(C{j})/2 -2*pi/numel(C{j})*(l-1) + 2*pi);
            Ccos{j}{l} = cos(mod(theta,pi)*2)*C{j}{l}; %TODO :Fix me
%             Ccos{j}{l}  = C{j}{l};
            Csin{j}{l} = sin(mod(theta,pi)*2)*C{j}{l};
        else
            Ccos{j}{l} = 0*C{j}{l};
            Csin{j}{l} = 0*C{j}{l};
        end
%         if j==1
%             C{j}{l} = 0*C{j}{l};
%         end
      end
    end
    
    disp('Inverse transform ...');
    %temp_restored = real(ifdct_wrapping(C,0,N,N));
    temp_restored = real(NifdctWrapping(C,is_real,M,N));
    temp_restored_cos = real(NifdctWrapping(Ccos,is_real,M,N));
    temp_restored_sin = real(NifdctWrapping(Csin,is_real,M,N));
    % L1 iterations
    
    for nupdate = 1:nell1
      disp(['Direct transform, within ell-1 iteration, shift nr.',num2str(n),' ...']);
      D = NfdctWrapping(temp_restored,is_real,finest,nbscales,nbangles_coarse);
      disp('Thresholding ...');
      thresh = nsigmas_coarse * sigma;
      for j = 1:length(C)
        if j == length(C), thresh = nsigmas_fine * sigma; end;
        for l = 1:length(C{j})
          thresh_jl = thresh*E{j}{l};
          
          % Uncomment for 'diagonal real thresholding'
          %                    redjl = real(D{j}{l});   
          %                    imdjl = imag(D{j}{l});
          %                    recjl = real(C{j}{l});
          %                    imcjl = imag(C{j}{l});
          %                    redjl = (recjl - redjl) .* (abs(recjl) > thresh_jl);
          %                    imdjl = (imcjl - imdjl) .* (abs(imcjl) > thresh_jl);
          %                    D{j}{l} = redjl + sqrt(-1)*imdjl;
          
          
          % Uncomment for 'diagonal complex thresholding'
          %                    modcjl = abs(C{j}{l});       
          %                    moddjl = abs(D{j}{l});          
          %                    argcjl = C{j}{l} ./ (modcjl+1e-16);
          %                    argdjl = D{j}{l} ./ (moddjl+1e-16);
          %                    moddjl = (modcjl - moddjl) .* (modcjl > thresh_jl);
          %                    D{j}{l} = argdjl .* moddjl;
          
          % Uncomment for 'block complex thresholding'
          modcjl = abs(C{j}{l});          
          moddjl = abs(D{j}{l});
          argcjl = C{j}{l} ./ (modcjl+1e-16);
          argdjl = D{j}{l} ./ (moddjl+1e-16);
          rowstep = M/N_rows{j}{l};
          %colstep = N/N_cols{j}{l};
          colstep = N/N_cols{j}{l};
          evenquad = ~mod(ceil(l*4/length(C{j})),2);
          if evenquad,
            if (j == 1)||(finest==2 && j==nbscales), fcolsjl = 1; else fcolsjl = F_cols{j}{l}; end;
            rowshift = - round(F_rows{j}{l}/fcolsjl * rowstep);
            testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[1 0]).^2 + circshift(modcjl,[-1 0]).^2 + ...
                                                      circshift(modcjl,[rowshift 1]).^2 + circshift(modcjl, [-rowshift -1]).^2));
          else
            if (j == 1)||(finest==2 && j==nbscales), frowsjl = 1; else frowsjl = F_rows{j}{l}; end;
            colshift = - round(F_cols{j}{l}/frowsjl * colstep);
            testcjl = sqrt(modcjl.^2 + neighb_weight*(circshift(modcjl,[0 1]).^2 + circshift(modcjl,[0 -1]).^2 + ...
                                                      circshift(modcjl,[1 colshift]).^2 + circshift(modcjl, [-1 -colshift]).^2));
          end;
          testcjl = testcjl ./ sqrt(1+4*neighb_weight*tuning_neighb);
          moddjl = (modcjl - moddjl) .* (testcjl > thresh_jl);
          D{j}{l} = argdjl .* moddjl;
          
        end
      end
      disp('Inverse transform ...')
      %temp_update = real(ifdct_wrapping(D,0,N,N));
      temp_update = real(NifdctWrapping(D,is_real,M,N));
      max(max(abs(temp_update)))
      
      temp_restored = temp_restored + temp_update;
    end;
    
    temp_restored = circshift(temp_restored,[-xshift, -yshift]);
    temp_restored_cos = circshift(temp_restored_cos,[-xshift, -yshift]);
    temp_restored_sin = circshift(temp_restored_sin,[-xshift, -yshift]);
    restored_img = (n-1)/n*restored_img + 1/n*temp_restored;
    
  end
end
restored_img = (restored_img(1:size(img,1),1:size(img,2)));

temp_restored_cos = temp_restored_cos(1:size(img,1),1:size(img,2));
temp_restored_sin = temp_restored_sin(1:size(img,1),1:size(img,2));
cos_old = temp_restored_cos;
sin_old = temp_restored_sin;
sum1 = sqrt(temp_restored_cos.^2 + temp_restored_sin.^2)+0.001;
temp_restored_cos = temp_restored_cos./sum1;
temp_restored_sin = temp_restored_sin./sum1;
angles = mod(atan2(sin_old,cos_old)+2*pi,2*pi)/2;
temp_restored_cos = cos(angles);
temp_restored_sin = sin(angles);


% restored_img = uint8(restored_img/max(restored_img(:))*255.0);
restored_img = uint8(restored_img);
% imwrite(restored_img,out_filename,'WriteMode','Append','Compression','None');
% imwrite(temp_restored_cos,out_filename_cos,'WriteMode','Append','Compression','None');
% imwrite(temp_restored_sin,out_filename_sin,'WriteMode','Append','Compression','None');
% 
%  [X,Y] = meshgrid(1:size(restored_img,1),1:size(restored_img,2));
% %  figure,imshow(restored_img);
% figure; quiver(Y',X',temp_restored_cos,temp_restored_sin,5);
%  figure,imagesc(temp_restored_cos);
%  figure, imagesc(temp_restored_sin);
% MSE = sum(sum((img-double(restored_img)).^2))/N^2;
% PSNR = 20*log10(255/sqrt(MSE));
% disp(['PSNR = ',num2str(PSNR)]);
% disp(['Time elapsed = ',num2str(toc)]);

% 
% figure(1),imshow(squeeze(restored_img*2));hold on;
% figure(1),imshow(squeeze(uint8(input/1.2)));hold on;
% [X,Y] = meshgrid(1:size(restored_img,2),1:size(restored_img,1));
% quiver(X,Y,double(restored_img)/50.*temp_restored_cos,-double(restored_img)/50.*temp_restored_sin);

% figure(2),imshow(squeeze(restored_img*2));hold on;
% figure(2),imshow(squeeze(input));hold on;
% [X,Y] = meshgrid(1:size(restored_img,2),1:size(restored_img,1));
% quiver(X,Y,double(restored_img)/50.*cos_old,-double(restored_img)/50.*sin_old);
% 
% figure(3),imagesc(temp_restored_cos);
% figure(4),imagesc(temp_restored_sin);


out_img  = restored_img;
dircos =  temp_restored_cos;
dirsin = temp_restored_sin;
gradmag = double(out_img);
Ix = 5*gradmag.*dirsin;
Iy = -5*gradmag.*dircos;
orientation = mod(atan2(dirsin,dircos),pi);
end



% figure(1); clf; imagesc(img); colormap gray; axis('image');
% title('Original Image')
% figure(2); clf; imagesc(noisy_img); colormap gray; axis('image');
% title('Noisy Image')
% figure(3); clf; imagesc(restored_img); colormap gray; axis('image');
% title('Restored image');
