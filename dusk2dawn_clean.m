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
% flds( cellfun(@(x) any(strcmp(x,{''})) , flds ) ) = []; % remove things that can't be varied  ? ONLY if I later make the line 29 bit more permissive for flexibility
for f = 1:length(flds)
    
    % get number of expected dimensions
    if     any( strcmp(flds{f}, {'asr_cutoff','asr_windowlength','ref_maxbadchannels', ...
            'chunk_len', 'ref_maxbadchannels', ...
            'chunk_overlap','asr_MaxMem','asr_maxdims','asr_UseRiemannian'}) ) % 1D variables
        expdim = 1;
    elseif any( strcmp(flds{f}, {'ref_tolerances' }) ) % 2D variables
        expdim = 2;
    else
        continue % skip any variables that can't be varied
%         expdim = 1; % assuming anything can be varied, but will need to be explicitly coded if expected dimension is not 1
    end               % ? note that this is dangerous with line 25 edit
    
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

%% handle varied parameters which affect whether to reuse the mask of relatively clean calibration data

if npars > 0 % if there are any varied pars then might need to reuse mask for reference data used for calibration
    ref_pars = cellfun( @(x)  startsWith(x,{'ref_','chunk_'}), pars.labels);  % if varied parameter affects choice of calibration data
    for p = 1:3 % loop through all 3 possible pars even if number of pars is less, because there is always a 3 par loop later
        if p <= npars 
            if ref_pars(p) 
                ref_dims(p) = size(pars.values{p},2); % if there is a varied par which affects ref data mask then have one element for each value 
            else
                ref_dims(p) = 1; % length of varied pars which don't determine the ref_mask is 1 because always the same mask
            end
        else
            ref_dims(p) = 1; % if there is no p'th varied par then just set that dim of ref_mask cell array to 1 anyway for easy indexing
        end
    end
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

    %% if computing multiple runthroughs with different parameters, then make empty cell array to store logical masks for reference data for this dataset
    if npars > 0       % so that ref data can be reused in multiple runthroughs
        ref_masks = cell(ref_dims(1),ref_dims(2),ref_dims(3));

        % if splitting by sleep stage then initialise the 2nd level cell arrays to store reused ref_masks for each stage
        if cfg.splitByStage % (if there are even multiple runthroughs we need to store calibration data for)
            ref_masks = cellfun( @(x) cell(1,nstages) ,ref_masks,'UniformOutput' ,false);
        end
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

                %% get indices for accessing reused calibration data masks
                if npars > 0
                    r1 = min(p1,size(ref_masks,1)); % if this varied parameter allows reusing the same ref_mask then index 1
                    r2 = min(p2,size(ref_masks,2)); % if not then index the correct element with the current param value 
                    r3 = min(p3,size(ref_masks,3)); % (min because if length of that ref_mask dimension is 1 then this will ensure index with 1)
                end

                %% check if splitting by stage and branch accordingly 
                if cfg.splitByStage

                    % loop through stages
                    for st = 1:nstages

                        % check for existing ref_mask for this set of pars & stage
                        if npars > 0
                            if ~isempty(ref_masks{r1,r2,r3}{st}) 
                                cfg2.ref_mask = ref_masks{r1,r2,r3}{st};
                            end
                        end
        
                        % run actual cleaning with the chosen parameters
                        [EEG_cleaned(st), cfg2_out(st)] = d2d_cleanData(EEG_stages(st), cfg2);

                        % store calibration mask in appropriate place (if not already done for this set of pars & stage)
                        if npars > 0 
                            if isempty(ref_masks{r1,r2,r3}{st}) 
                                ref_masks{r1,r2,r3}{st} = cfg2_out(st).ref_mask;
                            end
                        end
    
                        % run ICA for this stage
                        if cfg.ica.run
                            [EEG_cleaned(st), ~] = d2d_computeICA(EEG_cleaned(st), cfg);
                        end
                        
                        % run validation for this stage
                        [EEG_cleaned(st), ~] = d2d_validateCleaning(EEG_cleaned(st), cfg);
                        
                        fprintf('\n')
                    end
                    
                    % recombine stages of cleaned data
                    EEG_cleaned = d2d_recombineByStage(EEG_cleaned, EEG_dummy);
                    valid(count,:) = EEG_cleaned.etc.dusk2dawn.valid; % extract validation across stages for later merging

                    fprintf('\n\n')
                    
                else % not splitting by stages

                    % check for existing ref_mask for this set of pars
                    if npars > 0
                        if ~isempty(ref_masks{r1,r2,r3}) 
                            cfg2.ref_mask = ref_masks{r1,r2,r3};
                        end
                    end
    
                    % run actual cleaning with the chosen parameters
                    [EEG_cleaned, cfg2_out] = d2d_cleanData(EEG, cfg2);

                    % store calibration mask in appropriate place (if not already done for this set of pars)
                    if npars > 0 
                        if isempty(ref_masks{r1,r2,r3})
                            ref_masks{r1,r2,r3} = cfg2_out.ref_mask;
                        end
                    end
                    
                    % run ICA
                    if cfg.ica.run
                        [EEG_cleaned, ~] = d2d_computeICA(EEG_cleaned, cfg);
                    end
    
                    % run validation 
                    [EEG_cleaned, valid(count)] = d2d_validateCleaning(EEG_cleaned, cfg); 
    
                    fprintf('\n')
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
