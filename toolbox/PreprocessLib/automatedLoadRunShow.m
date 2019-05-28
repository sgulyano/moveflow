function [Iout,whatScale,whatScaleSums]=automatedLoadRunShow(depthImg)
tStart=tic;
type = ['a','b','s','n'];

for i=3
    tic
        
        options = struct('FrangiScaleRange',[1,10],'FrangiScaleRatio',1,'BlackWhite',false);
        [Iout,whatScale] = NFrangiFilter3D(depthImg,options,type(i));
        whatScaleSums=zeros(10,1);
    for j = 1:10
        whatScaleSums(j) = sum(sum(sum(whatScale==j)));
    end
        fname = strcat('Iout',type(i),'.mat')
        save(fname, 'Iout','whatScale','whatScaleSums');
        clear Iout whatScale whatScaleSums
    toc
end
toc(tStart);