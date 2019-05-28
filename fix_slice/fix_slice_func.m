function fix_slice_func(directory, filename, savedir, num_img, num_slice, switchTZ)
gfpinfo = imfinfo(fullfile(directory, sprintf(filename,0,0)));
H = gfpinfo.Height; W = gfpinfo.Width;
gfp = zeros(H, W, num_img+1, 'uint8');
pickslice = floor(num_slice/2);

hbar = waitbar(0, 'Please wait. This may take a while...');
% hbar = parfor_progressbar(num_img+1,'Please wait. This may take a while...');
for t = 0:num_img
    V = zeros(H, W, num_slice+1, 'uint8');
    Vori = zeros(H, W, num_slice+1, 'uint8');
    for z = 0:num_slice
        if switchTZ
            X = imread( fullfile(directory, sprintf(filename,z,t)) );
        else
            X = imread( fullfile(directory, sprintf(filename,t,z)) );
        end
        Vori(:,:,z+1) = X;
        % normalized
        X = single(X);
        X = X-min(X(:)); X = X/max(X(:)); X = uint8(round(X*255.0));
        
        % adjust intensity slice-wise
        n = 3;
        Idouble = im2double(X);
        avg = mean2(Idouble);
        sigma = std2(Idouble);
        X = imadjust(X,[max(avg-n*sigma,0) min(avg+n*sigma,1)],[]);
        
        V(:,:,z+1) = X;
    end
    % align slices
    tic;
    [Vnew, Valign] = align_stack(V, Vori, pickslice+1);
    toc;
    gfp(:,:,t+1) = uint8(max(Vnew,[],3)*255);
    
    for z = 0:num_slice
        if switchTZ
            imwrite( Valign(:,:,z+1), fullfile(savedir, sprintf(filename,z,t)) );
        else
            imwrite( Valign(:,:,z+1), fullfile(savedir, sprintf(filename,t,z)) );
        end
    end

    waitbar(t/num_img, hbar, ['Please wait. This may take a while... [' num2str(t) '/' num2str(num_img) ']']);
%     hbar.iterate(1); % update progress by one iteration 
end
imgflow2vdo( [savedir '_gfp.mp4'], gfp );
close(hbar); % close the progress bar
% uiwait(msgbox('Done'));
end