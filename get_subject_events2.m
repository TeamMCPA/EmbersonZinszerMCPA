function [event_matrix] = get_subject_events2(nirs_data,subject,channels,time_window,event_types)

%% Extract oxy data and marks from the Homer struct
oxy_timeser = nirs_data(subject).hmr.data.dConc(:,channels,1);
marks_vec = nirs_data(subject).otp.marksvector;

%% Grab windowed data according to onsets

% Build a matrix of all the onset marks. See temp_marks being revised to
% account for onset/offset marks. That line may need to be changed
% depending on the nature of the dataset being input.
marks_mat = nan(max(hist(marks_vec(marks_vec~=0))),length(event_types));
for type_i = 1:length(event_types),
    temp_marks = find(marks_vec==event_types(type_i));
    
    % Here is where we throw out the offsets
    temp_marks = temp_marks(1:2:end);
    
    marks_mat(1:length(temp_marks),type_i) = temp_marks;
end
% Throw out rows that are all NaNs
marks_mat = marks_mat(sum(isnan(marks_mat),2)<size(marks_mat,2),:);

%% Extract the individual events for each event type
% Set up matrices to store each channel's data for individual events
event_matrix = nan(length(time_window),length(channels),size(marks_mat,1),length(event_types));
    % event_matrix format:
    % event_matrix(time,channels,event_rep,event_type)

for type_i = 1:length(event_types),
    for event_j = 1:length(marks_mat(:,type_i)),
        if marks_mat(event_j,type_i) + time_window(end) <= length(oxy_timeser),
            event_matrix(:,:,event_j,type_i) = oxy_timeser(marks_mat(event_j,type_i)+time_window,:) - ones(length(time_window),1)*oxy_timeser(marks_mat(event_j,type_i)+time_window(1),:);
        else
            event_matrix(:,:,event_j:end,type_i) = NaN;
        end
    end
end
