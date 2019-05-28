function [ ddad_soma, ddad_dend, ddae_soma, ddae_dend ] = gui_select_neuron( V )
%GUI_SELECT_NEURON let user mark bounding box of dendrite and circle soma

num_neuron = 0;
ddad_soma = zeros(0,3);
ddad_dend = zeros(0,4);
ddae_soma = zeros(0,3);
ddae_dend = zeros(0,4);
hddad_box = [];
hddad_pnt = [];
hddae_box = [];
hddae_pnt = [];

fcn = makeConstrainToRectFcn('imrect',[1 size(V,2)], [1 size(V,1)]);

f = figure('Position',[100,100,1080,720]);

axes('Units','normalized','Position',[0.1,0.15,0.85,0.75],'FontSize',16);
imshow(V,[]);

hadd = uicontrol('Style','pushbutton','String','Add Neuron','FontSize',14,...
        'Units','normalized','Position',[0.15 0.05 0.3 0.05],...
        'Enable', 'on',...
        'Callback',@addCallback);

uicontrol('Style','pushbutton','String','Next Neuron','FontSize',14,...
        'Units','normalized','Position',[0.55 0.05 0.3 0.05],...
        'Enable', 'on',...
        'Callback',@nextCallback);

waitfor(f);
disp('Done');

    function addCallback(~,~)
        hadd.Enable = 'off';
        title('Draw ddaD Dendrite Bounding Box');
        hddad_box = imrect('PositionConstraintFcn', fcn);
        
        title('Circle ddaD Soma');
        hddad_pnt = imellipse;
        
        title('Draw ddaE Dendrite Bounding Box');
        hddae_box = imrect('PositionConstraintFcn', fcn);
        
        title('Circle ddaE Soma. Click Next to confirm selection.');
        hddae_pnt = imellipse;
    end

    function nextCallback(~,~)
        if ~isempty(hddad_box) && ~isempty(hddad_pnt) && ~isempty(hddae_box) && ~isempty(hddae_pnt) && ...
                isvalid(hddad_box) && isvalid(hddad_pnt) && isvalid(hddae_box) && isvalid(hddae_pnt)
            num_neuron = num_neuron + 1;
            ddad_bp = getPosition(hddad_box);
            ddad_dend(num_neuron,:) = [ddad_bp(1:2) ddad_bp(1:2)+ddad_bp(3:4)];
            delete(hddad_box);

            ddae_bp = getPosition(hddae_box);
            ddae_dend(num_neuron,:) = [ddae_bp(1:2) ddae_bp(1:2)+ddae_bp(3:4)];
            delete(hddae_box);

            cc = getPosition(hddad_pnt);
            cen = cc(1:2) + (cc(3:4)/2);
            rad = (cc(3)+cc(4))/4;
            ddad_soma(num_neuron,:) = [cen rad];
            delete(hddad_pnt);

            cc = getPosition(hddae_pnt);
            cen = cc(1:2) + (cc(3:4)/2);
            rad = (cc(3)+cc(4))/4;
            ddae_soma(num_neuron,:) = [cen rad];
            delete(hddae_pnt);

            hold on;
            plot(ddad_dend(num_neuron, [1 3 3 1 1]), ddad_dend(num_neuron, [2 2 4 4 2]), 'c--');
            viscircles(ddad_soma(num_neuron,1:2), ddad_soma(num_neuron,3), 'EdgeColor', 'g');
            plot(ddae_dend(num_neuron, [1 3 3 1 1]), ddae_dend(num_neuron, [2 2 4 4 2]), 'c--');
            viscircles(ddae_soma(num_neuron,1:2), ddae_soma(num_neuron,3), 'EdgeColor', 'm');
            hold off;

            title('Close figure when finish');
            hadd.Enable = 'on';
        else
            errordlg('Not both ddaD and ddaE are selected');
        end
    end
end

