function MCPA_struct = estimate_WindowAverages(data_file,incl_subjects,incl_channels,time_window)

%% estimate_WindowAverages takes HomER structs (with fields 'hmr' and 'opt') 
% for multiple subjects and prepares multichannel patterns for each 
% subject (one pattern for each condition discovered in the HomER file).
%
% The function is called with the following parameters:
% estimate_WidoowAverages( data_file , incl_subjects , incl_channels , time_window )
%
% data_file: a .mat file in HomER format containing all subjects'
% information in a size-n struct array for n subjects
%
% incl_subjects: a vector of indices for subjects to include in the
% analysis. Importantly the subject numbers correspond to the index in the
% struct array (e.g., MyData([1 3 5]) not any other subject number
% assignment.
%
% incl_channels: a vector of indices for channels to include in the
% analysis. Again, only the channel's position in the HomER struct matters,
% not any other channel number assignment.
%
% time_window: defined in number of measures. For data measured at 10 Hz,
% the time_window will be in 1/10 s units (e.g., 0:100 is 0-10 s). For data
% measured at 2 Hz, time will be in 1/2 s units (e.g., 0:20 is 0-10 s).
%
% The function will return a new struct containing some metadata and the 
% multichannel patterns for each participant and condition.

%% Load in the specified data file as "nirs_data" and then rename the data struct as "nirs_data"
fprintf('\nLoading data from %s...\n',char(data_file))
try
    nirs_data = load(data_file);
    data_name = fieldnames(nirs_data);
    fprintf('Found %.0f structs in the file:\n',length(data_name));
    fprintf('%s\n',data_name{:});
    fprintf('\nLoading first struct (%s)...',data_name{1});
    nirs_data = eval(sprintf('nirs_data.%s',data_name{1}));
    fprintf(' Done.\n');
catch
    fprintf(' Failed.\n');
    return
end

%% Scan through the marksvector and note all unique event types (marks)
fprintf('\nScanning for event types...',data_file)
event_types = [];
for subject = 1:length(incl_subjects)
    event_types = union(event_types,unique(nirs_data(incl_subjects(subject)).otp.marksvector));
end
% drop the initial zero type
event_types = event_types(2:end);
fprintf('\n%.0f event types found (numbered 1 to %.0f)\n',length(event_types),length(event_types));

%% Initialize the subj_mat matrix that will be output later (as NaNs for now)
% output format: subj_mat(event_type,channel,subject)
subj_mat = nan(...
    length(event_types),...
    length(incl_channels),...
    length(incl_subjects));

fprintf('\nA matrix size %.0f x %.0f x %.0f has been initialized.\nDimensions: (event types X channels X subjects)\n',size(subj_mat));


%% Extract data from each subject
fprintf('\nExtracting data for subject:\n');
for subject = 1:length(incl_subjects),
    fprintf(' %.0f',incl_subjects(subject));
    
    event_matrix = get_subject_events2(nirs_data,incl_subjects(subject),incl_channels,time_window,event_types);
    % event_matrix format:
    % event_matrix(time,channels,event_rep,event_type)
    
    event_means = nanmean(event_matrix,3);
    % event_means format:
    % event_means(time,channels,1,event_type)
    
    event_vects = nanmean(event_means,1);
    % event_vects format:
    % event_vects(1,channels,1,event_type)
    event_vects = reshape(event_vects,length(incl_channels),length(event_types));
    % event_vects format:
    % event_vects(channels,event_type)
    
    subj_mat(:,:,subject) = event_vects';
    % output format: subj_mat(event_type,channel,subject)
end
fprintf(' Done.\n');

fprintf('\nWriting MCPA_struct for this dataset...');
try
    MCPA_struct.data_file = data_file;
    MCPA_struct.created = datestr(now);
    MCPA_struct.summary_method = 'WindowAverages';
    MCPA_struct.time_window = time_window;
    MCPA_struct.incl_subjects = incl_subjects;
    MCPA_struct.incl_channels = incl_channels;
    MCPA_struct.patterns = subj_mat;
    fprintf(' Done.\n');
catch
    fprintf(' Failed.\n');
end
