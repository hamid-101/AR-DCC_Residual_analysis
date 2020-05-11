function [] = Step1_Modeling_AR_DCC(data)
addpath(genpath('Tools'));

current_dir = pwd;
DataFolder = [current_dir,'\Result'];

AR_ModelResultFolder = [DataFolder,'\1_whitendatabyAR\'];
GARCH_ModelResultFolder = [DataFolder,'\2_whitendatabyAR_GARCH\'];
DCC_ModelResultFolder = [DataFolder,'\3_AR_DCC\'];
       
    
for subj=1:length(data)
    clear data_ZeroConditionalMean

    %% Find Best AR and AR Modelling For removing Data Autocorrelation 
    mkdir(AR_ModelResultFolder);
    AR_ModelResult = [AR_ModelResultFolder,num2str(subj),'.mat'];
    if exist(AR_ModelResult,'file')
        load(AR_ModelResult)
    else        
         for i=1:size(data{subj},2)
            [dat2(:,i),~,~] = findbest_AR(data{subj}(:,i),5,0);%                
         end
         for i=1:size(data{subj},2)
             for j=i+1:size(data{subj},2)
                data_ZeroConditionalMean{i,j} = [dat2(:,i),dat2(:,j)];
             end
            fprintf('%d,',i);
        end
        save(AR_ModelResult,'data_ZeroConditionalMean');
    end
    %% Find Best GARCH Parameter
    BestArch=zeros(size(data{subj},2),size(data{subj},2),2);
    BestGarch=BestArch;
    data3 = cell(size(data{subj},2),size(data{subj},2));
    
    mkdir(GARCH_ModelResultFolder);
    GARCH_ModelResult = [GARCH_ModelResultFolder,num2str(subj),'.mat'];
    if exist(GARCH_ModelResult,'file')
        load(GARCH_ModelResult)
    else
        fprintf('Subject:%d\n',subj);
        clear dat3
        for i=1:size(data{subj},2)-1
            [dat3(:,i),BestArc(i),BestGarc(i)]=findbestGARCH_BIC(data_ZeroConditionalMean{i,i+1}(:,1),2);
            fprintf('%d,',i);
        end
        [dat3(:,i+1),BestArc(i+1),BestGarc(i+1)]=findbestGARCH_BIC(data_ZeroConditionalMean{i,i+1}(:,2),2);
            
        for i=1:size(data{subj},2)   
            for j=i+1:size(data{subj},2)
                data3{i,j} = [dat3(:,i),dat3(:,j)];
                BestArch(i,j,:)= [BestArc(i),BestArc(j)];
                BestGarch(i,j,:)= [BestGarc(i),BestGarc(j)];
            end
        end
        save(GARCH_ModelResult,'data3','BestArch','BestGarch');
    end
    
    %% Find Best DCC Parameter and DCC Modelling
    mkdir(DCC_ModelResultFolder);
    DCC_ModelResult = [DCC_ModelResultFolder,num2str(subj),'.mat'];
    if exist(DCC_ModelResult,'file')
        load(DCC_ModelResult)
    else
        result = cell(size(data{subj},2),size(data{subj},2));
        for i=1:size(data{subj},2)
            parfor j=i+1:size(data{subj},2)
            resAR = data_ZeroConditionalMean{i,j};   
            
            [result{i,j}.DCCParametr.dccP,...
                result{i,j}.DCCParametr.dccQ,...
                result{i,j}.DCCParametr.archP1,...
                result{i,j}.DCCParametr.archP2,...
                result{i,j}.DCCParametr.garchP1,...
                result{i,j}.DCCParametr.garchP2]=findbestDCC_BIC(resAR,BestArch(i,j,:),BestGarch(i,j,:));

            archP=[result{i,j}.DCCParametr.archP1 result{i,j}.DCCParametr.archP2];
            garchQ=[result{i,j}.DCCParametr.garchP1 result{i,j}.DCCParametr.garchP2];
            
            [parameters, ~, result{i,j}.Hthat, ~,  stdresid,...
            ~, ~, ~,~, ~,result{i,j}.Rthat]=...
            mydcc_mvgarch(resAR,result{i,j}.DCCParametr.dccP,result{i,j}.DCCParametr.dccQ,archP,garchQ);
            
            result{i,j}.DCCEstParametr=parameters;
            result{i,j}.DCCRes = stdresid;
           
        end
        fprintf('%d,',i);
       
        end
        save(DCC_ModelResult,'result');
    end
end
return;

function [resGARCH,BestArch,BestGarch]=findbestGARCH_BIC(y,maxlag)
if (nargin < 2)
    maxlag = 2;    
end
% tic
BestArch=0;
BestGarch=0;

% maxlag = 1;
ARCHEffect = zeros(maxlag+1,maxlag+1);
bicMatrix = ones(maxlag+1,maxlag+1)*inf;
aicMatrix = bicMatrix;
for archP=1:maxlag
    for garchP=1:maxlag
    try
        if archP==0 && garchP==0
            [H_rejection1, siglevel]=kstest(y);    
             ARCHEffect(archP+1 ,garchP+1)=siglevel;%[statistic;siglevel]';
             if (ARCHEffect(archP+1 ,garchP+1)>0.05)
                resGARCH = y;
                return;
             end
        else
            [parameters, likelihood, ~, ~, ht, ~] = fattailed_garch(y , archP ,garchP  , 'NORMAL');
                 [aicMatrix(archP+1 ,garchP+1),bicMatrix(archP+1 ,garchP+1)] = aicbic(likelihood,sum(abs(parameters)~=0),size(y,1));
            [~, siglevel]=kstest(y./sqrt(ht));    
             ARCHEffect(archP+1 ,garchP+1)=siglevel;%[statistic;siglevel]';
        end                            
        catch
            continue;
        end
    end
end
% toc
[a,b]= min(bicMatrix);
    [c,bestj] = min(a);
    besti = b(bestj)-1;
    bestj = bestj-1;
    BestArch=besti;
    BestGarch=bestj;
   
if ~(besti==0 && bestj==0)
    [parameters, likelihood, stderrors, robustSE, ht, scores] = fattailed_garch(y , besti ,bestj  , 'NORMAL');
     resGARCH=y./sqrt(ht);               
end


function [dccP,dccQ,archP1,archP2,garchP1,garchP2,archN,garchN]=findbestDCC_BIC(y,BestArch,BestGarch)
% ver 4
% correct error in temp<length
if (nargin <2) %% not good state,optimize all param
    maxlag =1;
    aic1 = ones(maxlag+1,maxlag+1,maxlag+1,maxlag+1,maxlag+1,maxlag+1)*inf;
    bic1 = ones(maxlag+1,maxlag+1,maxlag+1,maxlag+1,maxlag+1,maxlag+1)*inf;
    d = zeros(maxlag+1,maxlag+1,maxlag+1,maxlag+1,maxlag+1,maxlag+1);

    for dccP=0:maxlag
        for dccQ=0:maxlag
            for archP1=0:maxlag
                for archP2=0:maxlag
                    for garchP1=0:maxlag
                        for garchP2=0:maxlag
                            try
    if dccP==0 && dccQ==0
        continue;
    else
        archP=[archP1 archP2];%archP=zeros(1,size(y,2))
        garchQ=[garchP1 garchP2];%garchQ=zeros(1,size(y,2))
        [~, loglikelihood, ~, ~,  ~,...
        ~, ~, ~,~, ~,~]=mydcc_mvgarch(y,dccP,dccQ,archP,garchQ);
        % 
        [aic1(dccP+1,dccQ+1,archP1+1,archP2+1,garchP1+1,garchP2+1),bic1(dccP+1,dccQ+1,archP1+1,archP2+1,garchP1+1,garchP2+1)] =...
            aicbic(loglikelihood,sum([dccP,dccQ,archP,garchQ]),size(y,1));
        d(dccP+1,dccQ+1,archP1+1,archP2+1,garchP1+1,garchP2+1) = str2double(num2str([dccP,dccQ,archP1,archP2,garchP1,garchP2]')');
    end
                            catch
                                    continue;
                                end
                            end
                        end
                    end
                end
            end
    end
    b = reshape(d,1,[]);
    a = reshape(bic1,1,[]);
    [~,j] = min(a);

    temp= str2num(num2str(b(j))')';
    if length(temp)<length(size(bic1))
        temp = [0,temp];
    end
    dccP=temp(1);dccQ=temp(2);archP1=temp(3);archP2=temp(4);garchP1=temp(5);garchP2=temp(6);
else
    maxlag =1;
    aic1 = ones(maxlag+1,maxlag+1)*inf;
    bic1 = ones(maxlag+1,maxlag+1)*inf;
    d = zeros(maxlag+1,maxlag+1);

    for dccP=0:maxlag
        for dccQ=0:maxlag
            try
                if dccP==0 && dccQ==0
                    continue;
                else
                [~, loglikelihood, ~, ~,  ~,...
                ~, ~, ~,~, ~,~]=mydcc_mvgarch(y,dccP,dccQ,BestArch,BestGarch);
                % 
                [aic1(dccP+1,dccQ+1),bic1(dccP+1,dccQ+1)] =...
                    aicbic(loglikelihood,sum([dccP,dccQ])+1,size(y,1));
                d(dccP+1,dccQ+1) = str2double(num2str([dccP,dccQ]')');
            end
            catch
                    continue;
            end
        end
    end 
    b = reshape(d,1,[]);
    a = reshape(bic1,1,[]);
    [~,j] = min(a);

    temp= str2num(num2str(b(j))')';
    if length(temp)<length(size(bic1))
        temp = [0,temp];
    end
    dccP=temp(1);dccQ=temp(2);archP1=BestArch(1);archP2=BestArch(2);garchP1=BestGarch(1);garchP2=BestGarch(2);
    archN = BestArch;
    garchN = BestGarch;
end

function [resAR,BestAR,BestMA]=findbest_AR(y,maxlag,minlag)
resAR = y;
if (nargin < 3)
    minlag = 0;
end
if (nargin < 2)
    maxlag = 3;    
end
aic1 = ones(maxlag+1,maxlag+1)*inf;
bic1 = ones(maxlag+1,maxlag+1)*inf;
tic
    for j=minlag:maxlag
        for i=minlag:maxlag
            try
                if i==0 && j==0
                    results = lmtest1(y,5);
                    if sum(results.pval>0.05)==5
                        aic1(i+1,j+1)=-inf;
                        bic1(i+1,j+1)=-inf;
                    end                    
                    continue;
                else
                    [logL]=testarparam(y,i,j);
                    [aic1(i+1,j+1),bic1(i+1,j+1)] = aicbic(logL,i+j+2,size(y,1));
                end                            
            catch
                continue;
            end
        end
    end
    toc
    [a,b]= min(bic1);
    [c,bestj] = min(a);
    besti = b(bestj)-1;
    bestj = bestj-1;
    BestAR=besti;
    BestMA=bestj;
    
if ~(besti==0 && bestj==0)
    [~,resAR]=testarparam(y,besti,bestj);               
end

function [logL,res]=testarparam(y,i,j)
if i==0 && j~=0
    return;
    Mdl = arima('MALags',1:j);
elseif i~=0 && j==0
    Mdl = arima('ARLags',1:i);%1:i);
else
    return;
    Mdl = arima('ARLags',1:i,'MALags',1:j);%,'SARLags',1:k,'SMALags',1:m
    
end
[EstMdl,EstParamCov,logL1,info] = estimate(Mdl,y);
[res,~,logL] = infer(EstMdl,y);
   
