function [] = Show_DFC_Result()
current_dir = pwd;
DataFolder = [current_dir,'\Result'];

DFC_File = [DataFolder,'\6_DFC_AR_DCC\Result.mat'];
load(DFC_File);

Nstate = length(StateRate);

figure,
colormap jet
for StateIndex = 1:Nstate
    subplot(1,Nstate,StateIndex)
    ConMat = EstimatedDFCs{StateIndex};
    ConDur = StateRate(StateIndex);        
    imagesc(ConMat)
    axis square; axis ij 
    set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
    c = get(gca, 'Children');
    set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');        
%         text(1.5,-2,sprintf('Err:%0.2f Occure:%0.2f(s)',ConErr,ConDur/20), 'Fontsize', 10);
    title(sprintf('Time:%0.0f(s)',ConDur/(10*1/2)), 'Fontsize', 12);
end
