%% Dataset #1: A/V Discrimination
% data previously published in Emberson, Richards & Aslin (2015, control group)

% data is preprocessed in Homer2 such that we can save a *.nirs file for each subject after the preprocessing
% then each *.nirs file is loaded into MATLAB (along with timing
% information) and then loaded into a dataframe
% for more information on this procedure see https://sites.google.com/site/nirsanalysis/
data_file = 'alldata_WDIT_control2.mat';
incl_subjects = [1 3:6 8:14 17:18 20:21, 23:25];
incl_channels = [3:5,13:19];
time_window = [0:100];

condition1 = 1;
condition2 = 2;

% Use average dConc Oxy in a time window as point estimates for each
% channel for each event types and each subject (WDIT_window_average)
WDIT_window_averages = estimate_WindowAverages(data_file,incl_subjects,incl_channels,time_window);

% infant-level decoding
WDIT_infantleveldecoding = leave_one_Ss_out_classifyAverages(WDIT_window_averages,condition1,condition2);

% trial-level coding
WDIT_trialleveldecoding = leave_one_Ss_out_classifyIndividualEvents(WDIT_window_averages,condition1,condition2);

% examining infant-level decoding across different subset sizes for
% channels
setsize = 5;  
WDIT_infantleveldecoding_setsize5 = leave_one_Ss_out_classifyAverages(WDIT_window_averages,condition1,condition2,setsize);
WDIT_trialleveldecoding_setsize5 = leave_one_Ss_out_classifyIndividualEvents(WDIT_window_averages,condition1,condition2,setsize);

% now restrict the analyses to just the top 3 channels that were most
% informative and we have almost perfect performance for the 
incl_channels = [3:5,13:19];
incl_channels = incl_channels([1,3,8]);
WDIT_window_averages = estimate_WindowAverages(data_file,incl_subjects,incl_channels,time_window);
WDIT_infantleveldecoding_top2 = leave_one_Ss_out_classifyAverages(WDIT_window_averages,condition1,condition2);
WDIT_trialleveldecoding_top3 = leave_one_Ss_out_classifyIndividualEvents(WDIT_window_averages,condition1,condition2);