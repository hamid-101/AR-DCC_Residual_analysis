function []=Step3_ResidualAnalysis_AR_DCC(data)
current_dir = pwd;
DataFolder = [current_dir,'\Result'];

AR_DCC_Model_Folder = [DataFolder,'\3_AR_DCC\'];
Residual_Analysis_Folder = [DataFolder,'\5_Resid\'];
mkdir(Residual_Analysis_Folder);

NSubj = dir([AR_DCC_Model_Folder,'*.mat']);
for subj=1:length(NSubj) 
    fprintf('\n\tSubject : %d\n',subj);
    AR_DCC_Model = [AR_DCC_Model_Folder,num2str(subj),'.mat'];
    load(AR_DCC_Model);
    
    MyARDCC_result = result;
    
    
    Pvalue_be_Normal = zeros(size(MyARDCC_result,1),size(MyARDCC_result,2));
    Pvalue_be_UnCorrelated = Pvalue_be_Normal;
    Pvalue_be_Normal_InRaw = Pvalue_be_Normal;
    Pvalue_be_UnCorrelated_InRaw = Pvalue_be_Normal;
    ValidValues = Pvalue_be_Normal;
    
    for k=1:size(MyARDCC_result,1)
        for m=k+1:size(MyARDCC_result,2)
            testtype = 'spearman';            
            [Pvalue_be_Normal_InRaw(k,m),...
                Pvalue_be_UnCorrelated_InRaw(k,m)]=...
                Test_Normal_Uncorrelated(data{subj}(:,[k,m]),testtype);
            
            [Pvalue_be_Normal(k,m),...
                Pvalue_be_UnCorrelated(k,m)]=...
                Test_Normal_Uncorrelated(MyARDCC_result{k,m}.DCCRes,testtype);
            
            
            ValidValues(k,m) = 1;
        end
    end
    Residual_Analysis = [Residual_Analysis_Folder,num2str(subj),'.mat'];
    save(Residual_Analysis,'Pvalue_be_Normal','Pvalue_be_UnCorrelated',...
        'Pvalue_be_Normal_InRaw','Pvalue_be_UnCorrelated_InRaw',...
        'ValidValues');
end
return;

function [Pvalue_be_Normal,Pvalue_be_UnCorrelated]=Test_Normal_Uncorrelated(x,testtype)

%%% INDEPENDENCE among U residuals
M=size(x,1); % n. of signals
N=size(x,2); % length of signals
[~,Pvalue_be_UnCorrelated] = corr(x,'type',testtype);
Pvalue_be_UnCorrelated = Pvalue_be_UnCorrelated(1,2);
%% vectors for Jarque-Bera analysis
%% JARQUE-BERA TEST FOR GAUSSIANITY [from Luetkepohl 2005, New introduction to multiple time series analysis, pag 175]
U = x';
T=size(U,2); K=size(U,1); % notazione di Lutkepohl
Um=mean(U,2); % media di U
Su=zeros(K,K); % covarianza di U
for t=1:T
    Su=Su+(U(:,t)-Um)*(U(:,t)-Um)';
end
Su=Su./(T-1);
Ps=chol(Su)';
V = zeros(K,T);
for t=1:T
    V(:,t)=Ps\(U(:,t)-Um);
end
b1=zeros(K,1); b2=b1;
for k=1:K
    b1(k)=sum(V(k,:).^3) / T;
    b2(k)=sum(V(k,:).^4) / T;
end
%% statistica di Jarque-Bera (estensione a multivariate series)
lambdas=T*(b1'*b1)/6;
lambdak=T*(b2-3*ones(K,1))'*(b2-3*ones(K,1))/24;
dg=K;

Pvalue_be_Normal=1-chi2cdf(lambdas+lambdak,2*dg); %statistica totale (joint test)