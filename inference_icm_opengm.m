[version, executable, isloaded] = pyversion;
if ~strcmp(executable, '/Users/Yo/anaconda2/bin/python') && ~isloaded
    pyversion /Users/Yo/anaconda2/bin/python
end

%---------------------------------------------------------------
% MinSum  with ICM
%---------------------------------------------------------------
n=10;
nl=10;
unaries=py.numpy.random.rand(n , n, nl);
potts=py.opengm.PottsFunction(uint8([nl,nl]),0.0,0.05);
gm=py.opengm.grid2d2Order(unaries,potts);
%---------------------------------------------------------------
% Minimize
%---------------------------------------------------------------
%get an instance of the optimizer / inference-algorithm

inf=py.opengm.inference.Icm(gm);
% start inference (in this case verbose infernce)
inf.infer();
% get the result states
argmin=inf.arg();
% print the argmin (on the grid)
argmin.reshape(n,n)

result = double( py.array.array( 'd', py.numpy.nditer(argmin.reshape(n,n)) ) );
result = reshape(result, [n n])'; 