function [] = Step4_DFCAnalysis_AR_DCC()

current_dir = pwd;
DataFolder = [current_dir,'\Result'];

Correlation_Folder = [DataFolder,'\4_Corr\'];
DFC_Folder = [DataFolder,'\6_DFC_AR_DCC\'];
mkdir(DFC_Folder);

AR_DCCCorr = [];
NSubj = dir([Correlation_Folder,'*.mat']);
for subj=1:length(NSubj)
    Correlation_File = [Correlation_Folder,num2str(subj),'.mat'];
    load(Correlation_File); 
    AR_DCCCorr = cat(4,AR_DCCCorr,Rthat);
end

Ncluster = 0;
[EstimatedDFCs,StateRate,IndexofEachCorrelationtemp]=ClusterDynConnectivity(AR_DCCCorr,Ncluster);
    
save([DFC_Folder,'Result.mat'],'EstimatedDFCs','StateRate');


function [DFC_Cluster,rate,IDX,DFC_VarCluster]=ClusterDynConnectivity(DynamicCorrettions,Ncluster)
% Estimate Centroid by normal Distribution
% automatic Clustering when Ncluster is 0
% DynamicCorrettions n_ROI*n_ROI*Ntime*Nsubjects
MaxAutomaticElbowSearch = 12;
CUTOFFElbowSearch = 0.90;
[n_ROI a Ntime Nsubjects] = size(DynamicCorrettions);
for i=1:Nsubjects
    SelectTime{i} = (1:Ntime);
end
temp = triu(ones(n_ROI,n_ROI),1);
featureExtraction=@(X) abs((X(temp==1)));
Feature_Space = zeros(((n_ROI*(n_ROI-1)/2)),Ntime,Nsubjects);
DataClusterSpace = reshape(DynamicCorrettions,n_ROI,n_ROI,[]);
DynamicCorrettions(abs(DynamicCorrettions)<0.122)=0;
DataDistanceSpace = DynamicCorrettions;

for subj= 1:Nsubjects
    BestTime = SelectTime{subj};
    for i=1:length(BestTime)
        Feature_Space(:,i,subj) = featureExtraction(DataDistanceSpace(:,:,BestTime(i),subj));
    end
end
Feature_Space = reshape(Feature_Space,size(Feature_Space,1),[]);
Feature_Space(:,sum((Feature_Space))==-inf)=[];    
if Ncluster==0%% automatic clustering
     Feature_Space = zscore(Feature_Space,0,2);
%      tic, [IDX,C,sumd,Ncluster] = kmeans_opt(double(Feature_Space)',MaxAutomaticElbowSearch,CUTOFFElbowSearch,12);toc
    
    parfor i=2:MaxAutomaticElbowSearch
        [idx{i-1},~,sumd{i-1},D{i-1}]=kmeans(double(Feature_Space)',i,'Replicates',10,'Start',randi(3,[i,size(Feature_Space,1),10])-2);
        disp(i);
    end
    save('elbowkmeas1','idx','sumd','D');
    hold on
    within = zeros(MaxAutomaticElbowSearch-1,1);
    bithin = zeros(MaxAutomaticElbowSearch-1,1);
    for i=1:MaxAutomaticElbowSearch-1
    for n=1:size(idx{i},1)
        myindex = idx{i};
        myDist = D{i};
        within(i) = within(i)+myDist(n,myindex(n));
        Clusters = (1:size(myDist,2));
        Clusters(myindex(n))=[];
        bithin(i) = bithin(i)+sum(myDist(n,Clusters));
    end
    end
    figure(43),
    plot((2:MaxAutomaticElbowSearch),within./bithin,'-r');
    prompt = 'What is the Best Cluster? ';
    x = input(prompt);
    IDX = idx{x-1};
    Ncluster = x;
else
    Feature_Space = zscore(Feature_Space,0,2);
    [IDX,C,sumd,D] = kmeans(double(Feature_Space)',Ncluster,'Replicates',20,'Start',randi(3,[Ncluster,size(Feature_Space,1),20])-2);%,'Distance','correlation');
end
nIDX = zeros(Ntime,Nsubjects);
index = 1;
for subj= 1:Nsubjects
    BestTime = SelectTime{subj};
    nIDX(BestTime,subj) = IDX(index:index+length(BestTime)-1);
    index = index+length(BestTime);
end 
nIDX = reshape(nIDX,1,[]);
IDX = nIDX;
rate = zeros(Ncluster,1);
DFC_Cluster = cell(Ncluster,1);
DFC_VarCluster = cell(Ncluster,1);
for i=1:Ncluster
    Est_Rt = ((mean(DataClusterSpace(:,:,IDX==i),3)));
    Std_Rt = std(double(DataClusterSpace(:,:,IDX==i)),0,3);
    DFC_VarCluster{i} = Std_Rt;
    DFC_Cluster{i} = Est_Rt;
    rate(i)=sum(IDX==i);
end
[rate,I] = sort(rate,'descend');
DFC_Cluster = DFC_Cluster(I);
for i=1:Ncluster
IDX(IDX==I(i)) = i*10000;
end
IDX = IDX./10000;

DFC_VarCluster = DFC_VarCluster(I);
IDX = reshape(IDX,Ntime,Nsubjects);