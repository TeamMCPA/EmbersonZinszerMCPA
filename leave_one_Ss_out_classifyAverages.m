function allsubj_results = leave_one_Ss_out_classifyAverages(MCPA_pattern,cond1,cond2,setsize)  

%% leave_one_Ss_out_classifyAverages takes a struct output from the 
% estimate_WindowAverages function and performs the leave-one-subject-out
% test for each subject to classify the condition-averaged multichannel
% patterns.
%
% The function is called with the following parameters:
% leave_one_Ss_out_classifyAverages( MCPA_pattern , cond1 , cond2 [, setsize])
%
% MCPA_pattern: the output of the estimate_WindowAverages function. This
% structure should contain the metadata about how multichannel patterns
% were estimated (which subjects and channels were included) and the
% multichannel patterns themselves. 
% For backwards compatibility, a single 3D matrix may also be used here
% with dimensions: (event_types X channel X subject)
%
% cond1 and cond2: the event types to be compared (corresponding to the
% rows in MCPA_pattern.pattern, not necessarily the original mark number
% from the HomER file). If there are only 2 event types, cond1 and cond2
% may simply be 1 and 2, respectively. If there are 4 event types (e.g.,
% cat, dog, shoe, boot), then single conditions maybe be selected as a 
% subset (e.g., cond1=3 and cond2=4 for shoe vs. boot) or as vectors (e.g.,
% cond1=[1,2] and cond2=[3,4] for animals vs. footwear).
%
% setsize (optional): the number of channels that will be included in the
% analysis. If setsize is less than the total number of channels in
% MCPA_pattern struct (MCPA_pattern.incl_channels), then all combinations
% of the channels will be analized in a n-choose-k fashion. In this case,
% average decoding accuracy across combinations will be computed for each
% channel. Default behavior is to include all channels.


%% Preparing data for analysis
% Older version of WindowAverages only produced subj_mat instead of full
% MCPA_pattern struct, so analysis can be performed with only subj_mat data
% available, extracted from MCPA_pattern if necessary.

if isstruct(MCPA_pattern),
    
    if size(MCPA_pattern.patterns,2)~=length(MCPA_pattern.incl_channels),
        fprintf('\nWARNING: Size of patterns matrix is not consistent with number of channels in incl_channels.\nIgnoring incl_channels and proceeding with data in patterns matrix.\n\n');
    end
    subj_mat = MCPA_pattern.patterns;
    incl_subjects = MCPA_pattern.incl_subjects;
    incl_channels = MCPA_pattern.incl_channels;
    
elseif isa(MCPA_pattern,'float'),
    
    fprintf('\nWARNING: MCPA_pattern is class float.\nProceeding with floating point matrix data for backwards compatibility.\nParameters must be entered manually.\n\n');
    subj_mat = MCPA_pattern;
    incl_subjects = input('Please enter the vector of subjects to include (ex: [1 2 5]): ');
    incl_channels = input('Please enter the vector of channels to include (ex: [1 2 5]): ');
    data_file = input('Please enter path to the NIRS data file (*.mat): ');
    time_window = input('Please enter the time window (in scans) used to generate the multivoxel pattern (ex: [0:100]): ');
    
else
    
    fprintf('\nERROR: MCPA_pattern is not recognized as class struct or float.\n\n');
    return;
    
end


%% Check the setsize, if it is specified at all.

if nargin < 4
    % default for setsize is all channels
    setsize = size(subj_mat,2);
    
else
    
    if setsize > size(subj_mat,2)
        fprintf('\n ERROR: set size is larger than the total number of channels...\n\n');
        return;
    end
    
    fprintf('\n Running analysis based on a set size of %d channels...\n\n', setsize);

end



%% Report the test type and which conditions (event_types) are being compared

fprintf('\nLabel by Condition (2x2 similarity matrix for conditions [');
fprintf(' %d', cond1);
fprintf(' ] vs. [');
fprintf(' %d', cond2);
fprintf(' ])\n');


%% Setting up the combinations of channel subsets
% Create all possible subsets. If setsize is equal to the total number of
% channels, there will only be one 'subset' which is the full channel
% array. If setsize is less than the total number of channels, there will
% be n-choose-k subsets to analyze.

sets = nchoosek(1:size(subj_mat,2),setsize);
num_subjs = size(subj_mat,3);


%% Set up the results structure which includes a copy of MCPA_pattern

allsubj_results = [];

% save some basic info about in the allsubj_results
allsubj_results.MCPA_pattern = MCPA_pattern;
try
    allsubj_results.data_file = MCPA_pattern.data_file;
catch
    allsubj_results.data_file = data_file;
end
allsubj_results.created = datestr(now);
allsubj_results.test_type = 'Leave one subject out, Classify averaged conditions';
allsubj_results.setsize = setsize;
allsubj_results.incl_subjects = incl_subjects;
allsubj_results.incl_channels = incl_channels;
allsubj_results.cond1 = cond1;
allsubj_results.cond2 = cond2;
allsubj_results.subsets = incl_channels(sets);

% Initialize some nan matrices to store the results of all the
% combinations of channels when setsize < incl_channels
allsubj_results.accuracy.subjXchan = nan(num_subjs,size(subj_mat,2));
allsubj_results.accuracy.subsetXsubj = nan(size(sets,1),num_subjs);



%% Iterate through all the subjects and attempt classification

if setsize==size(subj_mat,2), 
    fprintf('  Subject\tAcc\tSimilarity Matrix\n');
else
    fprintf('  Subject-wise results suppressed for subset analysis.\n');
    %fprintf('  Channel\tMean Acc\n');
end

for i = 1:num_subjs,

    % Group holds the _indexes_ of the included subjects, not their actual
    % subject numbers. Drop the index for the subject who is going to be
    % tested in this iteration of the loop (i is his/her index)
    group = [1:num_subjs];
    group(i) = [];
    
    % Build the averaged data for the group (non-test) subjects
    group_vecs = mean(subj_mat(:,:,group),3)';
    group_vecs = cat(2,mean(group_vecs(:,cond1),2),mean(group_vecs(:,cond2),2));
    
    % Extract the averaged data for the test subject
    test_vecs = subj_mat(:,:,i)';
    test_vecs = cat(2,mean(test_vecs(:,cond1),2),mean(test_vecs(:,cond2),2));
    
    % To store the results for the number of channels, by the number of
    % channels, by whether or not the channel was accurate for that search
    % light as well as the accurcies for the matching vs. the mismatching
    % conditions
    results = nan(size(sets,1),size(subj_mat,2));
    
    for curr_subset = 1:size(sets,1)
        
        %% Pull out the group model and test vectors and compare them
        if setsize >= 2          
            % For most setsize values, we can use correlation of the
            % channel-wise vectors to test similarity
            
            g1t1 = corr(group_vecs(sets(curr_subset,:),1),test_vecs(sets(curr_subset,:),1));
            g2t2 = corr(group_vecs(sets(curr_subset,:),2),test_vecs(sets(curr_subset,:),2));
            sort1 = g1t1 + g2t2;
            
            g1t2 = corr(group_vecs(sets(curr_subset,:),1),test_vecs(sets(curr_subset,:),2));
            g2t1 = corr(group_vecs(sets(curr_subset,:),2),test_vecs(sets(curr_subset,:),1));
            sort2 = g1t2 + g2t1;
            
        elseif setsize == 1
            % If setsize is equal to one, then we just take the differences
            % between the model and test values. This is not a very ideal
            % approach and shouldn't be used, but it's here to catch the
            % situations where somebody enters setsize 1
            
            g1t1 = abs(group_vecs(sets(curr_subset,:),1)-test_vecs(sets(curr_subset,:),1));
            g2t2 = abs(group_vecs(sets(curr_subset,:),2)-test_vecs(sets(curr_subset,:),2));
            sort1 = mean([g1t1,g2t2]);
            
            g1t2 = abs(group_vecs(sets(curr_subset,:),1)-test_vecs(sets(curr_subset,:),2));
            g2t1 = abs(group_vecs(sets(curr_subset,:),2)-test_vecs(sets(curr_subset,:),1));
            sort2 = mean([g1t2,g2t1]);   
        end
        accuracy = sort1>sort2;
        
        
        %% Report and store results
        if setsize == size(subj_mat,2) 
            % When setsize is equal to the full channel array (only one
            % subset to test), the results can be printed directly to the
            % Matlab command window
            if accuracy,
                fprintf('  %.0f\t\t%d\t%0.2f + %0.2f > %0.2f + %0.2f\n',incl_subjects(i),accuracy,g1t1,g2t2,g1t2,g2t1)
            else
                fprintf('  %.0f\t\t%d\t%0.2f + %0.2f < %0.2f + %0.2f\n',incl_subjects(i),accuracy,g1t1,g2t2,g1t2,g2t1)
            end
            results(curr_subset,sets(curr_subset,:)) = accuracy;
            
        elseif setsize == 1
            results(curr_subset,sets(curr_subset,:)) = accuracy;
            
        else
            results(curr_subset,sets(curr_subset,:)) = accuracy;

        end
        
            
    end
    
    % now save it in a struct across subjects
    foo =(squeeze(nanmean(results,1)));
    allsubj_results.accuracy.subjXchan(i,:) = foo(:);
    setresults.acc = nan(size(sets,1),1);
    
    % save the result for each set (we only look at the first channel
    % because the decoding accuracy is the same for each channel in the
    % set)
    for j = 1:size(sets,1)
        setresults.acc(j) = results(j,sets(j,1));

    end
    allsubj_results.accuracy.subsetXsubj(:,i) = setresults.acc;
    
    allsubj_results.accuracy.subsetXchanXsubj(:,:,i) = results(:,:);
    
end


%% If subsets were used, display figures

if setsize ~= size(subj_mat,2) 
    
    % Channel-wise mean accuracy figure
    figure;
    if isstruct(MCPA_pattern), 
        xvals = MCPA_pattern.incl_channels;
    else
        xvals = 1:size(allsubj_results.accuracy.subjXchan,2);
    end
    errorbar(1:length(xvals),mean(allsubj_results.accuracy.subjXchan),std(allsubj_results.accuracy.subjXchan)/sqrt(size(allsubj_results.accuracy.subjXchan,1)));
    title('Decoding Accuracy across all channels')
    set(gca,'XTick',[1:length(xvals)])
    set(gca,'XTickLabel',xvals)
    
    % Subject-wise mean accuracy figure
    figure
    errorbar(1:length(incl_subjects),mean(allsubj_results.accuracy.subjXchan'),repmat(std(mean(allsubj_results.accuracy.subjXchan'))/sqrt(size(allsubj_results.accuracy.subjXchan,2)),1,size(allsubj_results.accuracy.subjXchan,1)));
    title('Decoding Accuracy across all subjects')
    set(gca,'XTick',[1:length(incl_subjects)])
    set(gca,'XTickLabel',incl_subjects)
    
end



