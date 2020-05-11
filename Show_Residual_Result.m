function [] = Show_Residual_Result()
current_dir = pwd;
DataFolder = [current_dir,'\Result'];

Residual_Analysis_Folder = [DataFolder,'\5_Resid\'];
AllPvalue_be_Normal = [];
AllPvalue_be_UnCorrelated = [];

AllPvalue_be_Normal_InRaw = [];
AllPvalue_be_UnCorrelated_InRaw = [];

NSubj = dir([Residual_Analysis_Folder,'*.mat']);
for subj=1:length(NSubj) 
    fprintf('\n\tSubject : %d\n',subj);
    Residual_Analysis = [Residual_Analysis_Folder,num2str(subj),'.mat'];
    load(Residual_Analysis);
    AllPvalue_be_Normal = cat(3,AllPvalue_be_Normal,Pvalue_be_Normal(ValidValues==1));
    AllPvalue_be_UnCorrelated = cat(3,AllPvalue_be_UnCorrelated,Pvalue_be_UnCorrelated(ValidValues==1));
    
    AllPvalue_be_Normal_InRaw = cat(3,AllPvalue_be_Normal_InRaw,Pvalue_be_Normal_InRaw(ValidValues==1));
    AllPvalue_be_UnCorrelated_InRaw = cat(3,AllPvalue_be_UnCorrelated_InRaw,Pvalue_be_UnCorrelated_InRaw(ValidValues==1));
end
figure,
X1 = reshape(AllPvalue_be_Normal,1,[]);
X2 = reshape(AllPvalue_be_Normal_InRaw,1,[]);
X = [X1;X2]';
[count,center] = hist(X,20);
XGuassian = [center,count./repmat(sum(count),size(count,1),1)];
bar(XGuassian(:,1),XGuassian(:,[2,3]));
legend('AR-DCC','Raw Data');
xlabel('P-values');
ylabel('frequency of the probability ');
set(gca,'FontName','Times New Roman','FontSize',12)
title('Histogram of P values in the normality test');
ylim([0 0.9])

figure,
X1 = reshape(AllPvalue_be_UnCorrelated,1,[]);
X2 = reshape(AllPvalue_be_UnCorrelated_InRaw,1,[]);
X = [X1;X2]';
[count,center] = hist(X,20);
XUncorrelated = [center,count./repmat(sum(count),size(count,1),1)];
bar(XUncorrelated(:,1),XUncorrelated(:,[2,3]));
xlabel('P-values');
ylabel('frequency of the probability ');
set(gca,'FontName','Times New Roman','FontSize',12)
title('Histogram of P values in the non-correlated test');
ylim([0 0.9])

