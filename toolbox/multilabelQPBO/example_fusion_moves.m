addpath('qpboMex-master');

%% example 1: binary
nNodes=4;
%source,sink
terminalWeights=[
    0,16;
    0,13;
    20,0;
    4,0
];

%From,To,Capacity,Rev_Capacity
edgeWeights = zeros(2,2,5);
edgeWeights(:,:,1) = [0 10; 4 0];
edgeWeights(:,:,2) = [0 12; -1 0];
edgeWeights(:,:,3) = [0 -1; 9 0];
edgeWeights(:,:,4) = [0 14; 0 0];
edgeWeights(:,:,5) = [0 0; 7 0];

edge_ends=[
    1,2;
    1,3;
    2,3;
    2,4;
    3,4
    ];

[cut, labels] = fusion_moves(terminalWeights, edgeWeights, edge_ends, 5);
fprintf('Example 1: Cut = %d, labels = [%d; %d; %d; %d]\n', cut, labels);
% correct answer: cut = 22; labels = [1; 1; 2; 1];
if ~isequal(cut, 22)
    warning('Wrong value of cut!')
end
if ~isequal(labels, [1; 1; 2; 1])
    warning('Wrong energy of labels!')
end

%% example 2: multi-label
nNodes=4;
%source,sink
terminalWeights=[
    10,0,16;
    0,7,13;
    20,11,0;
    4,0,1
];

%From,To,Capacity,Rev_Capacity
edgeWeights = zeros(3,3,5);
edgeWeights(:,:,1) = [0 0 0; 0 0 10; 4 0 0];
edgeWeights(:,:,2) = [0 0 0; 0 5 12; -1 0 0];
edgeWeights(:,:,3) = [0 0 -1; 0 0 0; 9 0 0];
edgeWeights(:,:,4) = [0 14 0; 0 0 0; 0 0 0];
edgeWeights(:,:,5) = [0 0 0; 0 0 0; 7 0 0];

edge_ends=[
    1,2;
    1,3;
    2,3;
    2,4;
    3,4
    ];

[cut, labels] = fusion_moves(terminalWeights, edgeWeights, edge_ends, 5);
fprintf('Example 2: Cut = %d, labels = [%d; %d; %d; %d]\n', cut, labels);
% correct answer: cut = 10; labels = [1; 1; 3; 3];
if ~isequal(cut, 10)
    warning('Wrong value of cut!')
end
if ~isequal(labels, [1; 1; 3; 3])
    warning('Wrong energy of labels!')
end 