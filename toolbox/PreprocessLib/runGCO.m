function [minCut,Eall,Dall,Sall,Lall]=runGCO(depthImg,tau,type)

intDepthImg=int32(depthImg);
dDepthImg=double(depthImg);
[rows,cols,depth]=size(depthImg);
minCut=zeros(size(depthImg));

aParam=.8;
if nargin<2
    tau=20;
    type='2D';
end
if strcmp(type,'2D')
    pixels = rows*cols;
    for page=1:depth
        disp(num2str(page));
        img=intDepthImg(:,:,page);
        h=GCO_Create(pixels,2);
        dataCost=zeros(rows,cols,2);
        weights=zeros(size(dataCost));
        dataCost(:,:,1)=img;
        dataCost(:,:,2)=2*tau-img;
        dataCostRow1=reshape(dataCost(:,:,1),1,pixels);
        dataCostRow2=reshape(dataCost(:,:,2),1,pixels);
        dataCost=cat(1,dataCostRow1,dataCostRow2);
%         dataCost=aParam.*dataCost;
        dataCost=int32(dataCost);
        GCO_SetDataCost(h,dataCost);
        clear dataCost

        img=dDepthImg(:,:,page);
        weights(:,:,1)=img-circshift(img,[0,-1]);

        weights(:,:,2)=img-circshift(img,[-1,0]);
        weights=exp(-abs(weights));
        weights(:,cols,1)=0;%there should be no associated weight between the last edge and the first edge
        weights(rows,:,2)=0; %No associated weight from the bottom edge to the top
        weights=int32(weights);
        weights=double(weights);

        iIndices=cat(2,1:pixels,1:pixels);
        jIndices=cat(2,(1+rows):pixels,1:rows,2:pixels,1); %Set up the indices for the sparse matrix construction; much faster than iterating through for loops
        values=reshape(weights,pixels,2);
        clear weights;
        values=values';
        nborWeights=sparse(iIndices,jIndices,values,pixels,pixels);
        clear values;

        GCO_SetNeighbors(h,nborWeights);
        clear nborWeights;
        GCO_Expansion(h);                % Compute optimal labeling via alpha-expansion 
        labeledList = GCO_GetLabeling(h);
%         [E D S L]=GCO_ComputeEnergy(h);
%         Eall(page)=E;
%         Dall(page)=D;
%         Sall(page)=S;
%         Lall(page)=L;

        binaryImg=reshape(labeledList,rows,cols);
        minCut(:,:,page)=binaryImg;
        clear labeledList binaryImg
        clear gco_matlab
    end
    minCut=minCut-1; %change labels 1 and 2 to 0 and 1
elseif strcmp(type,'3D')
voxels = numel(depthImg);
    h=GCO_Create(voxels,2);
    dataCost=zeros(rows,cols,depth,2);
    weights=zeros(size(dataCost));
    dataCost(:,:,:,1)=intDepthImg;
    dataCost(:,:,:,2)=2*tau-intDepthImg;
    dataCostRow1=reshape(dataCost(:,:,:,1),1,voxels);
    dataCostRow2=reshape(dataCost(:,:,:,2),1,voxels);
    dataCost=cat(1,dataCostRow1,dataCostRow2);
%     dataCost=aParam.*dataCost;
    dataCost=int32(dataCost);
    GCO_SetDataCost(h,dataCost);
    clear dataCost
    
    
    weights(:,:,:,1)=dDepthImg-circshift(dDepthImg,[0,-1,0]);
    weights(:,:,:,2)=dDepthImg-circshift(dDepthImg,[-1,0,0]);
    weights(:,:,:,3)=dDepthImg-circshift(dDepthImg,[0,0,-1]);
    weights=exp(-abs(weights));
    weights(:,cols,:,1)=0;%there should be no associated weight between the last edge and the first edge
    weights(rows,:,:,2)=0; %No associated weight from the bottom edge to the top
    weights(:,:,depth,3)=0;%No associated weight from the front page to the back
    weights=int32(weights); %must be an integer for the C++ code, but
    weights=double(weights);%must be double for sparse matrix construction
    
    iIndices=cat(2,1:voxels,1:voxels,1:voxels);
    jIndices=cat(2,(1+rows):voxels,1:rows,2:voxels,1,(1+(rows*cols)):voxels,1:rows*cols); %Set up the indices for the sparse matrix construction; much faster than iterating through for loops
    values=reshape(weights,voxels,3);
    clear weights;
    values=values';
    nborWeights=sparse(iIndices,jIndices,values,voxels,voxels);
    clear values;

    
    
%     for i=1:rows
%         for j=1:cols
%             x=(j-1)*rows+i;
%             y=x+rows;
%             nborWeights(x,y)=weights(i,j,1);
%             y=x+1;
%             nborWeights(x,y)=weights(i,j,2);
%         end
%     end
    GCO_SetNeighbors(h,nborWeights);
    clear nborWeights;
    GCO_Expansion(h);                % Compute optimal labeling via alpha-expansion 
    labeledList = GCO_GetLabeling(h);
    [E D S L]=GCO_ComputeEnergy(h);
    Eall=E;
    Dall=D;
    Sall=S;
    Lall=L;
    
    minCut=reshape(labeledList,rows,cols,depth);
    minCut=minCut-1; %change labels 1 and 2 to 0 and 1
    clear labeledList
    clear gco_matlab
else
    disp('Unknown GCO option.  Valid options are ''2D'' or ''3D''');
end
