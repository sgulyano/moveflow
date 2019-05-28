function output=runEpicflow(im1,im2,matchfn,setting,t)
	path=pwd;
    if nargin == 5,
        sedfile = ['debug/edge1_' num2str(t) '.bin'];
        im1file = ['debug/im1_' num2str(t) '.png'];
        im2file = ['debug/im2_' num2str(t) '.png'];
        outfile = ['debug/output_' num2str(t) '.flo'];
    else
        sedfile = 'debug/edge1.bin';
        im1file = 'debug/im1.png';
        im2file = 'debug/im2.png';
        outfile = 'debug/output.flo';
    end
    runfile = ['./external/epicflow/epicflow-static ' im1file ' ' im2file ' ' sedfile ' ' matchfn ' ' outfile ' -' setting];
    
	SED(im2uint8(im1), fullfile(path,sedfile));
	cd(path);
	imwrite(im1,im1file);
	imwrite(im2,im2file);
	system(runfile);
	output=readFlowFile(outfile);
end