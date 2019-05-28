[version, executable, isloaded] = pyversion;
if ~strcmp(executable, '/Users/Yo/anaconda2/bin/python') && ~isloaded
    pyversion /Users/Yo/anaconda2/bin/python
end
%%
img = imread('cameraman.tif');

shape = size(img);%[12,12];
numVar = shape(1)*shape(2);

data = im2double(img);%round(rand(shape), 1);
labelsA = double(data>0.5);

beta = 0.5;
gm = py.opengm.gm( 2*ones(1,numVar,'int32') );

% add unaries
disp('Add unaries...');
st = tic;
unaries = ones([2 shape]);
unaries(1,:,:)=data;
unaries(2,:,:)=1.0-data;

unaryFunctionIds = gm.addFunctions(py.numpy.array(unaries(:)').reshape(-1,2));
gm.addFactors(unaryFunctionIds,py.numpy.arange(numVar));
stats.unary = toc(st);

disp('Add Potts regularizer...');
st = tic;
pottsFunction = py.opengm.pottsFunction(int16([2,2]),0.0,beta);
pottsFunctionId = gm.addFunction(pottsFunction);

for x = 0:shape(1)-1
   for y = 0:shape(2)-1
        if x+1 < shape(1)
            vi0 = y+x*shape(2);
            vi1 = y+(x+1)*shape(2);
            gm.addFactor(pottsFunctionId,int32([vi0,vi1]));
        end
        if y+1 < shape(2)
            vi0 = y +x*shape(2);
            vi1 = y+1 + x*shape(2);
            gm.addFactor(pottsFunctionId,int32([vi0,vi1]));
        end
   end
end
stats.potts = toc(st);

disp('Add Field of Expert...');
st = tic;
block4Function = zeros([2,2,2,2]);
block4Function(1,1,1,1)=2.0;
%block4Function(2,2,2,2)=10.0;
block4FunctionId = gm.addFunction(py.numpy.array(block4Function(:)').reshape(2,2,2,2));

for x = 0:shape(1)-1
   for y = 0:shape(2)-1
        if x+1 < shape(1) && y+1 < shape(2)
            vi0 = y + x*shape(2);
            vi1 = y+1 + x*shape(2);
            vi2 = y + (x+1)*shape(2);
            vi3 = y+1 + (x+1)*shape(2);
            vis = [vi0,vi1,vi2,vi3];
            gm.addFactor(block4FunctionId,int32(vis));
        end
   end
end
stats.foe = toc(st);

%%
disp('Inference using BP');
st = tic;
parameter = py.opengm.InfParam(pyargs('steps',1000,'damping',0.9,'convergenceBound',0.001));
inf2 = py.opengm.inference.BeliefPropagation(gm,pyargs('parameter',parameter));

inf2.infer();
arg=inf2.arg();
labelsB = arg.reshape(shape);
stats.infer = toc(st);
% disp(labelsA)

labelsB_row = double( py.array.array('d', py.numpy.nditer(labelsB)) );  % Add order='F' to get data in column-major order (as in Fortran 'F' and Matlab
labelsB_size = cell2mat(cell(labelsB.shape));
labelsB_mat = reshape(labelsB_row, labelsB_size); 
% disp(labelsB_mat)

figure, imshow(labelsB_mat);