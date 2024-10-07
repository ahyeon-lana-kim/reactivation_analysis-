function getReps_ROI()

    % Set the correct path where the load_nii function is
    niiPath = '' 

    % Path for functional data - you'll need to create this folder and put all
    % clean MNI files there (the path name suits my PC, but may be different on
    % yours)
    func_path = '';

    % Get a list of functional files that end with .nii, if your files end with .gz, in the terminal run gunzip *.nii.gz
    files = dir(fullfile(func_path, '*.nii'));

    % get the file names into an array
    ss_list = strings(length(files), 1);
    for i = 1:length(files)
        ss_list(i) = files(i).name;
    end

    % The path where the ROI is (change if different)
    ROIpath = '';

    % Output path (create this folder)
    pathO = '';

    % Accounting for the HRF delay
    iLag = 0; % 
    hrfLag = 4; 

   % Load the boundary file
    load(fullfile(func_path, 'allBoundaries.mat'));

    % Load boundary file (create a variable containing the TRs of all event boundaries. Then save this variable as a allBoundaries.mat file, and you can load it again whenever you run this function.)
    sceneBegin = [1; allBoundaries(1:end-1) + 1]; % beginning of scenes
    DiffBoundaries = allBoundaries - sceneBegin; % lengths

    numEvents = length(allBoundaries); % the number of events in the movie

    % Set the path to the ROI directory
    ROIpath = '';

    % Create list of ROI names
    ROIfiles = dir(fullfile(ROIpath, 'ROI.nii'));


    % For cluster use - get the argument for this iteration (number of ROI)
    %jobID = '1';for debugging only. REMOVE this line when running on the cluster
    %iROI = str2num(jobID);

    iROI=1;

    % Get name of ROI
    ROIname = ROIs{iROI};

    % Load ROI
    nii = load_nii(fullfile(ROIpath, ROIname));

    % Extract ROI data
    ROI = nii.img;

    % Ensure that ROI data has only zeroes and ones
    %hist(ROI(:))

    % Generate the number of ROI voxels
    voxIndices = find(ROI == 1);
    numVox = length(voxIndices); % Correctly define numVox as the number of ROI voxels

    % Print the number of voxels for debugging
    disp(['numVox: ', num2str(numVox)]);

    % Ensure numEvents is defined and scalar
    if ~exist('numEvents', 'var') || ~isscalar(numEvents)
        error('numEvents is not defined or not a scalar.');
    end

    % Preallocate space for boundary representations
    EB = zeros(numVox, numEvents, length(files));
    scene = EB;
    CTRp = EB; % past controls
    CTRf = EB; % future controls

    % Initialize idx for ROI timecourses
    [x, y, z] = ind2sub(size(ROI), voxIndices);
    idx = [x, y, z];

   
        % Load the functional data
        %nii = load_nii(fullfile(func_path, files(iSub).name));
        %data = nii.img;
        % Define the path and filename

   % Define the path
    func_path = '';

    % Get a list of all .nii files in the directory
    files = dir(fullfile(func_path, '*.nii'));

    % Loop through each file
    for iFile = 1:length(files)
    % files(iFile).name refers to the current file's name in the iteration
    % Construct the full file path
    full_func_file_path = fullfile(func_path, files(iFile).name);

    % Load the functional data
    nii = load_nii(full_func_file_path);
    data = nii.img;

    % Now 'data' contains the image data from the NIfTI file
    disp(['Processing file: ', files(iFile).name]); % Display the current file name
    disp(size(data)); % Display the size of the data for verification
    
    % You can add further processing steps here for each 'data'
    end


        % Get ROI timecourses
        tcs = zeros(length(idx), size(data, 4));
        for j = 1:length(idx)
            tcs(j, :) = zscore(double(data(idx(j, 1), idx(j, 2), idx(j, 3), :)));
        end
        data = tcs'; clear tcs

        % Get representations, while correcting for the HRF lag
        for iEvent = 1:length(allBoundaries)
            if allBoundaries(iEvent) + iLag > size(data, 1) % if we can't correct, because there isn't enough data
                EB(:, iEvent, iSub) = nan;
                scene(:, iEvent, iSub) = nan;
                continue
            else
                eventBound = allBoundaries(iEvent);
                % Extract the corrected event boundary representations
                EB(:, iEvent, iSub) = data(eventBound + iLag + hrfLag, :);

                if eventBound + 10 < size(data, 1)
                    CTRf(:, iEvent, iSub) = data(eventBound + iLag + hrfLag + 10, :);
                end

                if eventBound - 10 > 0
                    CTRp(:, iEvent, iSub) = data(eventBound + iLag + hrfLag - 10, :);
                end

                % Extract the averaged scene representation. Remove 5 TRs from
                % the beginning and end of each scene, to separate scenes from
                % boundaries.
                scene(:, iEvent, iSub) = mean(data(sceneBegin(iEvent) + 5 + iLag + hrfLag:allBoundaries(iEvent) - 5 + iLag + hrfLag, :), 1);
            end
        end
        clear nii data

        disp(['Finished processing file ', num2str(iSub)]);
    end

    % Save representations
    EB(EB == 0) = nan;
    scene(scene == 0) = nan;
    CTRp(CTRp == 0) = nan;
    CTRf(CTRf == 0) = nan;

    save(fullfile(pathO, ['six_' ROIname(1:end-4) '_lag' num2str(iLag) '_reps']), 'EB', 'scene', 'CTRp', 'CTRf');

    disp('Processing completed and results saved.');
    
    end 


