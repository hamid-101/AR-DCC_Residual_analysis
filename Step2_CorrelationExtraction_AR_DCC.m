function [] = Step2_CorrelationExtraction_AR_DCC()

current_dir = pwd;
DataFolder = [current_dir,'\Result'];

AR_DCC_ModelResultFolder = [DataFolder,'\3_AR_DCC\'];
Correlation_Folder = [DataFolder,'\4_Corr\'];
mkdir(Correlation_Folder);

NSubj = dir([AR_DCC_ModelResultFolder,'*.mat']);
for subj=1:length(NSubj)
    AR_DCC_ModelResult = [AR_DCC_ModelResultFolder,num2str(subj),'.mat'];
    load(AR_DCC_ModelResult);
    
    Rthat = ones(size(result,1),size(result,1),size(result{1,2}.Rthat,3));
    for i=1:size(result,1)
        for j=i+1:size(result,1)
            Rthat(i,j,:) = (squeeze((result{i,j}.Rthat(1,2,:))));
            Rthat(j,i,:) = Rthat(i,j,:);
        end
    end
    
    
    Correlation_File = [Correlation_Folder,num2str(subj),'.mat'];
    save(Correlation_File,'Rthat');
end