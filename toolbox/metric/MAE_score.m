function score = MAE_score( goldswc, testswc )
%EV_MAE compute mean absolute error (MAE) of the obtained trace against the
%manually acquired ground truth. 
%MAE = mean_i(min_j|p_i-q_j|) + mean_i(min_k|q_i-p_k|)

if ischar(goldswc)
    goldswc = read_swc_file(goldswc);
end
if ischar(testswc)
    testswc = read_swc_file(testswc);
end

goldswc = resample_neuron(goldswc);
testswc = resample_neuron(testswc);

gmin = zeros(1,size(goldswc,1));
for i = 1:size(goldswc,1)
    gmin(i) = min(sqrt(sum(bsxfun(@minus, testswc(:,3:5), goldswc(i,3:5)).^2,2)));
end

tmin = zeros(1,size(testswc,1));
for i = 1:size(testswc,1)
    tmin(i) = min(sqrt(sum(bsxfun(@minus, goldswc(:,3:5), testswc(i,3:5)).^2,2)));
end

score = (mean(gmin) + mean(tmin)) / 2;
end

