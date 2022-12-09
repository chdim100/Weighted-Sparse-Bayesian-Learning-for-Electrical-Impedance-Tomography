%%%% Loading Open 2D EIT data %%%%
%
% Aku Seppänen
% University of Eastern Finland
% 17.3.2017
%

% If you wish to use only part of data (not all 79 current injections), select as follows (examples):
USE_ONLY_PART_PART_OF_DATA = 1;
InjectionSet = 4; % Change this to use different sets.

if InjectionSet == 1
    UseTheseInj = 1:16; Injname = 'adjacent';
elseif InjectionSet == 2
    UseTheseInj = 17:32; Injname = 'skip1';
elseif InjectionSet == 3
    UseTheseInj = 33:48; Injname = 'skip2';
elseif InjectionSet == 4
    UseTheseInj = 49:64; Injname = 'skip3';
elseif InjectionSet == 5
    UseTheseInj = 65:79; Injname = 'all_against1';
elseif InjectionSet == 6
    UseTheseInj = [17:32 49:64]; % combination of skip1 and skip3
end

% Specify Case:
savedatapath = 'data_mat_files';
datasetno = 2; objectno = 5;  % this is one example only (corresponds to Case2.5)
datamat_file = [savedatapath,'/datamat_',num2str(datasetno),'_',num2str(objectno)]

% Then load the data and pick parts of data accordingly:
load(datamat_file)
Ninj = size(CurrentPattern,2); % Total number of current injections
Nmeas = size(MeasPattern,2); % Number of measurements per current injection
if USE_ONLY_PART_PART_OF_DATA
   CurrentPattern = CurrentPattern(:,UseTheseInj);
   Uel = Uel(:,UseTheseInj);
   Ninj = length(UseTheseInj); % Number of used current injections
end
Uel = Uel(:);

