
load('SampleData\SimData.mat');
%% Processing
Step1_Modeling_AR_DCC(data);
Step2_CorrelationExtraction_AR_DCC();
Step3_ResidualAnalysis_AR_DCC(data);
Step4_DFCAnalysis_AR_DCC();
%% Show Result
%DFC
Show_DFC_Result();
%residual
Show_Residual_Result();
