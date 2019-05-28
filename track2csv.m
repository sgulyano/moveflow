function track2csv( savefile, soma_pos, num_timestep, a_sum, a_mu, a_num )
%TRACK2CSV save track results to CSV file
T = table((0:num_timestep)');
T.Properties.VariableNames{end} = 'ImageNumber';
Ts = cell(1,length(soma_pos));
for i = 1:length(soma_pos)
    Ts{i} = table(soma_pos{i}(:,1), soma_pos{i}(:,2));
    Ts{i}.Properties.VariableNames = {['Xsoma_' num2str(i)], ['Ysoma_' num2str(i)]};
end
Tsum = table(a_sum);
Tsum.Properties.VariableNames{1} = 'IntensitySum';
Tmu = table(a_mu);
Tmu.Properties.VariableNames{1} = 'IntensityMean';
Tnum = table(a_num);
Tnum.Properties.VariableNames{1} = 'PixelCount';
writetable([T Ts{:} Tsum Tmu Tnum], savefile);
disp(['Save track result to ' savefile ' ... Done']);
end

