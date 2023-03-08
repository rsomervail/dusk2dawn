%  
%   Usage: 
%     -> call this function via the GUI: pop_dusk2dawn_clean(EEG)
%  
%   Clean dataset with 1 or more sets of ASR parameters & compute validation metrics
%   (note that output EEG structure is the raw data, and cleaned data must be loaded 
%    with d2d_loadData)
% 
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function EEG = dusk2dawn_clean(EEG, cfg)

%% handle multiple files
nfiles = length(EEG);
EEG_all = EEG; % store all datasets as other variable, doing this even if only one dataset for simplicity

% get all filenames in group (even if only one)
allFiles = {EEG.filename};

%% handle varied parameters
flds = fields(cfg); pars.labels = {}; pars.values = {};
for f = 1:length(flds)
    
    % get number of expected dimensions
    if     any( strcmp(flds{f}, {'asr_cutoff','asr_windowlength','ref_maxbadchannels','chunk_len','chunk_overlap','asr_MaxMem'}) ) % 1D variables
        expdim = 1;
    elseif any( strcmp(flds{f}, {'ref_tolerances' }) ) % 2D variables
        expdim = 2;
    else
        continue % skip any variables that can't be varied
    end
    
    % check if variable was varied
    if sum( size( cfg.(flds{f}) )>1 ) == expdim
        pars.labels{end+1} = flds{f};
        if expdim == 1
            pars.values{end+1} =  mat2cell( cfg.(flds{f}) , size(cfg.(flds{f}),1) , ones(1,size(cfg.(flds{f}),2)) )  ;
        elseif expdim == 2
            cfg.(flds{f}) = cfg.(flds{f})';
            pars.values{end+1} =  mat2cell( cfg.(flds{f}) , size(cfg.(flds{f}),1) , ones(1,size(cfg.(flds{f}),2)) )  ;
        end
    elseif sum( size( cfg.(flds{f}) )>1 ) > expdim
        error('pop_dusk2dawn: too many dimensions of input for parameter ''%s''', flds{f})
    end
end
npars = length(pars.labels);
if npars > 3
    error('pop_dusk2dawn: cannot support more than 3 varied parameters: %s', strjoin( pars.labels,', ') )
end

% re-order varied parameters so that the raw data will be plotted next to the less severe ASR level
for p = 1:npars
    label = pars.labels{p};
    vals = pars.values{p};
    vals = [vals{:}]; % arrange into array to allow sorting
    if any( strcmp(label, {'asr_cutoff','ref_maxbadchannels'}) )
        [~,sortIndy] = sort( vals ,'descend'); 
    elseif strcmp(label, {'ref_tolerances'})  % assuming ref_tolerances are 2x1 cells
        [~,sortIndy] = sort( range(vals) ,'descend'); 
    else
        [~,sortIndy] = sort( vals ,'ascend'); 
    end
    pars.values{p} = pars.values{p}(sortIndy);
end

% store parameters for later
cfg.pars = pars; 
% get lengths of any varied parameters
if npars >= 1, np1 = size(pars.values{1},2); else, np1 = 1; end
if npars >= 2, np2 = size(pars.values{2},2); else, np2 = 1; end
if npars == 3, np3 = size(pars.values{3},2); else, np3 = 1; end

%% check whether it's appropriate to store mask containing data for calibration and reuse in each iteration
if npars > 0 % this only matters if there are multiple varied parameters
    % if not varying calibration data parameters & no sliding window
    if     ~any(startsWith(pars.labels, 'ref_')) && ~cfg.splitBySlidingWindow  
        flag_reuseCalibData = true;
    else
        flag_reuseCalibData = false;
    end
else
    flag_reuseCalibData = false;
end

%% generate all filenames that cleaned data will be saved as later ____________________________________________________________
count = 1;
% parameter loop 1
for p1 = 1:np1

    % get first varied parameter (if exists)
    if npars >= 1
%         saveas1 = [ '_' pars.labels{1} '_' num2str(p1) ]; % using parameter name
%         saveas1 = [ '_p1-' num2str(p1,'%02d') '-' num2str(np1,'%02d') ]; % using code
        saveas1 = [ '_p1-' num2str(p1,'%02d')  ]; % using simplest code
    else
        saveas1 = '';
    end

    % parameter loop 2
    for p2 = 1:np2

        % get second varied parameter (if exists)
        if npars >= 2
%             saveas2 = [ '_' pars.labels{2} '_' num2str(p2) ]; % using parameter name
%             saveas2 = [ '_p2-' num2str(p2,'%02d') '-' num2str(np2,'%02d') ]; % using code
            saveas2 = [ '_p2-' num2str(p2,'%02d')  ]; % using simplest code
        else
            saveas2 = '';
        end

        % parameter loop 3
        for p3 = 1:np3

            % get third varied parameter (if exists)
            if npars == 3
%                 saveas3 = [ '_' pars.labels{3} '_' num2str(p3) ];
%                 saveas3 = [ '_p3-' num2str(p3,'%02d') '-' num2str(np3,'%02d') ]; % using code
                saveas3 = [ '_p3-' num2str(p3,'%02d') ]; % using simplest code
            else
                saveas3 = '';
            end

            % store filename for later saving
            cleanfile_suffixes{count} = [ '_clean' saveas1 saveas2 saveas3 ];

            count = count + 1;
        end % p1
    end % p2
end % p3
%% _________________________________________________________________________________________________________________________

%% LOOP THROUGH DATASETS
tIN = tic;
for f = 1:nfiles

    %% get dataset
    EEG = EEG_all(f);
    cur_setname = EEG.setname;
    
    % if data not loaded already, then load the data & store this info so it can be cleared later to preserve memory
    if ischar(EEG.data)
        dataflag = EEG.data;
        temp = pop_loadset('filepath',EEG.filepath,'filename',EEG.filename);
%         temp = pop_loadset('eeg',EEG);
        EEG.data = temp.data; clear temp
    else
        dataflag = '';
    end

    fprintf([ repmat('_', 1,200) '\n\n']);
    fprintf([ '%s: cleaning & validating dataset: "%s" (%d/%d)\n'], mfilename, EEG.setname,f,nfiles)
    fprintf([ repmat('_', 1,200) '\n\n']);

    %% first check the data is appropriate for ASR-cleaning
    EEG = d2d_checkData(EEG);
    
    %% store all dusk2dawn settings and info in EEG structure
    cfg.origFile   = [ strrep( EEG.setname,'.set','') '.set']; % easier than conditionals
    cfg.cleanFiles = cellfun( @(x) [ strrep(cfg.origFile,'.set','') x '.set' ]  , cleanfile_suffixes , 'UniformOutput'  , false);
    EEG.etc.dusk2dawn.cfg = cfg;

    EEG.etc.dusk2dawn.group.datasets = allFiles; % store other files of this group
    
    % split data by stage
    if cfg.splitByStage
        [EEG_stages, EEG_dummy] = d2d_splitByStage(EEG, cfg); % output split data stages and also the original data 
        nstages = length(EEG_stages);                         % structure without the data (dummy) for later recombination with d2d_combineByStage
        if nstages > length(cfg.stageCodes)
           EEG_dummy.etc.dusk2dawn.cfg.stageCodes{end+1} = '[unscored]'; 
        end
        clear EEG % save RAM by clearing non-split version of the data 
    end
    
    %% run validation on raw data for reference
    fprintf([ repmat('-', 1,200) '\n\n']);
    fprintf('%s: computing validation metrics on raw data as a baseline\n\n',mfilename)
    fprintf([ repmat('-', 1,200) '\n\n']);

    % loop through sleep stages 
    if cfg.splitByStage
        for st = 1:nstages
    
            % record any varied parameters here for later validation
            if npars > 0
                for p = 1:npars
                    EEG_stages(st).etc.dusk2dawn.asr.pars.(['par_' pars.labels{p}]) = {nan};
                end
            end
    
            EEG_stages(st).etc.dusk2dawn.asr.raw = true;
    
            % run ICA
            if cfg.ica.run
                [EEG_stages(st), ~] = d2d_computeICA(EEG_stages(st), cfg);
            end
            
            % run validation for this dataset
            [EEG_stages(st), ~] = d2d_validateCleaning(EEG_stages(st), cfg);
    
        end
    
        % recombine raw data & save (even if not validating)
        EEG = d2d_recombineByStage(EEG_stages, EEG_dummy); 
        valid_orig = EEG.etc.dusk2dawn.valid;
    
        % save raw dataset
        ctemp = [];
        ctemp.saveName = cfg.origFile;
        ctemp.savePath = cfg.savePath;
        EEG = d2d_saveDataset(EEG, ctemp); 
        clear EEG % delete afterwards to save memory while dataset is split by stage
    
    else % not splitting by stages
    
        % record any varied parameters here for later validation
        if npars > 0
            for p = 1:npars
                EEG.etc.dusk2dawn.asr.pars.(['par_' pars.labels{p}]) = {nan};
            end
        end
    
        EEG.etc.dusk2dawn.asr.raw = true;
    
        % run ICA
        if cfg.ica.run
            [EEG, ~] = d2d_computeICA(EEG, cfg);
        end
    
        % run validation 
        [EEG, valid_orig] = d2d_validateCleaning(EEG, cfg); 
    
        % save raw dataset
        ctemp = [];
        ctemp.saveName = cfg.origFile;
        ctemp.savePath = cfg.savePath;
        EEG = d2d_saveDataset(EEG, ctemp); % keep if not splitting by stage
        
    end % end if split by stages

    fprintf(['\n' repmat('-', 1,200) '\n\n']);
    fprintf('%s: computed validation metrics on raw data as a baseline\n',mfilename)
    fprintf('total time elapsed = %.1f mins\n',toc(tIN)/60)
    fprintf(['\n' repmat('-', 1,200) '\n\n']);
    
    %% loop through all varied parameters & clean data
    cfg2 = rmfield( cfg, 'savePath' ); % make copy for ASR cleaning which shouldn't save the files until after validation
    clearvars cfg2_out
    count = 1;

    % loop 1
    for p1 = 1:np1
    
        % get first varied parameter (if exists)
        if npars >= 1
            cfg2.(pars.labels{1}) = pars.values{1}{p1}; % set the actual parameter
            cfg2.pars.(['par_' pars.labels{1}]) = {pars.values{1}{p1}};  % also store it for later validation
        end
        
        % loop 2
        for p2 = 1:np2
    
            % get second varied parameter (if exists)
            if npars >= 2
                cfg2.(pars.labels{2}) = pars.values{2}{p2}; % set the actual parameter
                cfg2.pars.(['par_' pars.labels{2}]) = {pars.values{2}{p2}};  % also store it for later validation
            end
            
            % loop 3
            for p3 = 1:np3
    
                % get third varied parameter (if exists)
                if npars == 3
                    cfg2.(pars.labels{3}) = pars.values{3}{p3}; % set the actual parameter
                    cfg2.pars.(['par_' pars.labels{3}]) = {pars.values{3}{p3}};  % also store it for later validation
                end

                % print 
                if npars > 0
                    fprintf(['\n' repmat('-', 1,200) '\n']);
                    fprintf([ '%s: cleaning dataset: "%s" - runthrough: %02d/%02d   -   START'  '\n'],mfilename, cur_setname, count, np1*np2*np3 );
                        val = pars.values{1}{p1};
                        if size(val,1) > 1, val = val'; end % if 2D param then transpose before converting to str
                        str =       sprintf('%s = %s (%02d/%02d)',      pars.labels{1}, num2str(val),  p1, np1);
                    if npars >= 2
                        val = pars.values{2}{p2};
                        if size(val,1) > 1, val = val'; end % if 2D param then transpose before converting to str
                        str = [ str sprintf(', %s = %s (%02d/%02d)',    pars.labels{2}, num2str(val),  p2, np2) ];
                    end
                    if npars == 3
                        val = pars.values{3}{p3};
                        if size(val,1) > 1, val = val'; end % if 2D param then transpose before converting to str
                        str = [ str sprintf(' and %s = %s (%02d/%02d)', pars.labels{3}, num2str(val),  p3, np3) ];
                    end
                    fprintf([ 'varied parameters: ' str ],mfilename)
                    fprintf(['\n' repmat('-', 1,200) '\n\n']);
                end
                
                % loop through sleep stages
                if cfg.splitByStage
                    for st = 1:nstages

                        % use already-generated calibration mask if it is compatible
                        if flag_reuseCalibData
                            if count > 1 % if not first file
                                cfg2.ref_mask = cfg2_out(st).ref_mask;
                            end
                        end
        
                        % run actual cleaning with the chosen parameters
                        [EEG_cleaned(st), cfg2_out(st)] = d2d_cleanData(EEG_stages(st), cfg2);
    
                        % run ICA for this stage
                        if cfg.ica.run
                            [EEG_cleaned(st), ~] = d2d_computeICA(EEG_cleaned(st), cfg);
                        end
                        
                        % run validation for this stage
                        [EEG_cleaned(st), ~] = d2d_validateCleaning(EEG_cleaned(st), cfg);
    
                    end
                    
                    % recombine stages of cleaned data
                    EEG_cleaned = d2d_recombineByStage(EEG_cleaned, EEG_dummy);
                    valid(count,:) = EEG_cleaned.etc.dusk2dawn.valid; % extract validation across stages for later merging
                    
                else % not splitting by stages

                    % use already-generated calibration mask if it is compatible
                    if flag_reuseCalibData
                        if count > 1 % if not first file
                            cfg2.ref_mask = cfg2_out.ref_mask;
                        end
                    end
    
                    % run actual cleaning with the chosen parameters
                    [EEG_cleaned, cfg2_out] = d2d_cleanData(EEG, cfg2);
                    
                    % run ICA
                    if cfg.ica.run
                        [EEG_cleaned, ~] = d2d_computeICA(EEG_cleaned, cfg);
                    end
    
                    % run validation 
                    [EEG_cleaned, valid(count)] = d2d_validateCleaning(EEG_cleaned, cfg); 
    
                end % end if split by stages
    
                % save cleaned dataset
                ctemp = [];
                ctemp.savePath = cfg.savePath;
                ctemp.saveName = cfg.cleanFiles{count}; 
                d2d_saveDataset(EEG_cleaned, ctemp);

                % print 
                if npars > 0
                    fprintf(['\n' repmat('-', 1,200) '\n']);
                    fprintf([ '%s: cleaning dataset: "%s" - runthrough: %02d/%02d   -   END'  '\n'],mfilename, cur_setname, count, np1*np2*np3 );
                    fprintf('total time elapsed = %.1f mins\n',toc(tIN)/60)
                    fprintf(['\n' repmat('-', 1,200) '\n\n']);
                end

                % increment counter
                count = count + 1;
                
            end % loop 3    
        end % loop 2
    end % loop 1
    
    %% merge results of post-ASR validation across datasets and store in all datasets
    ctemp = [];
    if cfg.splitByStage
        ctemp.valid = [ valid_orig; valid];
    else
        ctemp.valid = [ valid_orig, valid]';
    end
    ctemp.loadPath  = cfg.savePath;
    ctemp.loadNames = cfg.cleanFiles;
%     ctemp.loadNames = cellfun( @(x) [ x '.set' ], cfg.cleanFiles, 'UniformOutput', false);
    ctemp.loadNames = [cfg.origFile, ctemp.loadNames ];
    ctemp.pars = cfg.pars;
    valid_merged = d2d_validateMerge( ctemp  );
    
    % recombine raw data if it was split into stages
    if cfg.splitByStage
        EEG = d2d_recombineByStage(EEG_stages, EEG_dummy);
    end
    
    % output raw data with validation structure added  (no need to resave header with this info because its done in the d2d_validateMerge function)
    EEG.etc.dusk2dawn.valid_merged = valid_merged;

    % reinsert into EEG_all
    if ~isempty(dataflag)
        EEG.data = dataflag;
    end
    EEG_all(f) = EEG;

    % print dataset output
    fprintf(['\n\n' repmat('_', 1,200) '\n\n']);
    fprintf( '%s: finished cleaning & validating dataset: "%s" (%d/%d)\n', mfilename, EEG.setname,f,nfiles)
    fprintf('total time elapsed = %.1f mins\n',toc(tIN)/60)
    fprintf([ repmat('_', 1,200) '\n\n\n']);

end % loop through datasets

% output set of datasets appropriately processed
EEG = EEG_all; % rename to EEG so it can be passed back out to master function

% merge results of validation across the group
if nfiles > 1
    EEG = d2d_group_validateMerge(EEG);
end

fprintf('\n\n')
fprintf('%s: finished cleaning all selected datasets\n',mfilename)
fprintf('total time elapsed = %.1f mins\n',toc(tIN)/60)
fprintf(['\n\n' repmat('_ _', 1,100) '\n\n']);
fprintf('\n\n')

end % function
