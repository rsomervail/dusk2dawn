% 
%   - Merges the results of post-cleaning validation across runs and raw data
% 
% Dusk2Dawn
% Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022
%           www.iannettilab.net 
%%  
function [valid_merged] = d2d_validateMerge( cfg )
fprintf('merging validation metrics across ASR parameters ...\n')

%% get filenames and paths for each dataset whose validation structures will be merged
% if no datasets specified then make user select the datasets with validation structures to merge
if ~isfield(cfg,'loadNames') 
    [cfg.loadNames, cfg.loadPath] = uigetfile('*.set', 'choose file to load',  'MultiSelect', 'on');
end

%% get validation structure & list of varied parameters
if isfield(cfg,'valid') && isfield(cfg,'pars') % if varied pars & validation structures to merge have been provided
    valid = cfg.valid;
    pars = cfg.pars;
else % otherwise load both from file
        
    % load validated dataset (info/header only) 
    for f = 1 : length(cfg.loadNames)
        EEG = pop_loadset('filepath', cfg.loadPath, 'filename', cfg.loadNames{f}, 'loadmode', 'info');
        cfg.valid(f,:) = EEG.etc.dusk2dawn.valid; % ! check this works
    end 
    valid   = cfg.valid;
    pars    = EEG.etc.dusk2dawn.cfg.pars;
   
end
npars = length(pars.labels);

%% reorder validation structure so that original (uncleaned file) comes first (necessary for handling multiple dimensions later)
[~, sortInd] = sort( [valid(:,1).raw], 'descend' );
valid = valid(sortInd,:); clear sortInd

%% add raw data to each varied parameter
if npars > 0
    for p = 1:npars
        pars.values{p} = [ {nan}, pars.values{p}];
    end
end

%% merge 1st order metrics
% loop through stages 
nstages = size(valid,2);
for st = 1:nstages

    %% get sleep stage if split by stages
    if nstages > 1 
        stage = valid(1,st).stage;
    end

    %% loop through datasets  
    nfiles = size(valid,1);
    for f = 1 : nfiles
        
    %% load validation subfield for this file from EEG in memory or disk, or from cfg directly
    v = valid(f,st);
    
    %% for this file, get the index for each parameter value that was varied 
    if npars > 0
        parindy = nan(1,npars);
        for p = 1:npars
            temp = find(  cellfun( @(x) isequaln(x, v.asr.pars.(['par_' pars.labels{p}]){1} ), pars.values{p} )   );
            if ~isempty(temp)
                parindy(p) = temp;
            else
                error 'd2d_validateMerge: chosen parameters include a value that is not present in any selected dataset'
            end
        end % parindy is then a list of indexes for this file f, which tells later sections where to place the 
            % validation metric values in the output matrices
        pstr = sprintf('parindy(%d),',1:npars); % turn into string so it can be passed to eval throughout the code        
        pstr(end) = []; % remove trailing comma 
    else
        pstr = num2str(f);
    end
    
    %% general ASR stats 
    eval([ 'asr.m_propCleaned('         pstr ') = v.asr.propCleaned;' ])
    eval([ 'asr.m_propGFPRemoved('      pstr ') = v.asr.propGFPRemoved;' ])
    eval([ 'asr.m_propVarRemoved('      pstr ') = v.asr.propVarRemoved;' ])
    eval([ 'asr.m_tElapsed('            pstr ') = v.asr.tElapsed;' ])
    
    % calibration data stats
    temp  = v.asr.calibDataLen;
    eval([ 'asr.m_calibDataLen_min('    pstr ') = min(temp);' ])
    temp2 = 100 * sum(temp > 60) / length(temp);
    eval([ 'asr.m_calibDataLen_prop60s('    pstr ') = temp2;' ])
    clear temp temp2
    
    %% fft
    if isfield(v,'fft')
        if v.fft.run
            fft.binFreqs        =  v.fft.binFreqs;
            fft.binFreqsLabels  =  v.fft.binFreqsLabels;
            eval([ 'fft.m_binAmp('           pstr ',:,:)      = v.fft.binAmp;' ])   
            eval([ 'fft.m_binAmp_avgChan('   pstr   ',:)      = v.fft.binAmp_avgChan;' ])
        end
    end
    
    %% SW
    if isfield(v,'sw')
        if v.sw.run
            eval([ 'sw.m_nwaves('                pstr     ')   = v.sw.nwaves;'       ])
            eval([ 'sw.m_amp_mean('         pstr     ')   = v.sw.amp_mean;'     ])
            eval([ 'sw.m_amp_median('       pstr     ')   = v.sw.amp_median;'   ])
            eval([ 'sw.m_tval_median('      pstr     ')   = v.sw.tval_median;'  ])
    %         eval([ 'sw.waves('              pstr ',:,:)   = v.sw.waves;' ])
    %         eval([ 'sw.waves_median(' pstr   ',:)   = v.sw.waves_median;'  ])
                      
        end
    end
    
    %% ICA
    if isfield(v,'ica')
        if v.ica.run
        
        % general ICA stuff
        eval([ 'ica.m_numIC('    pstr ') = v.ica.numIC;' ])
        eval([ 'ica.m_tElapsed(' pstr ') = v.ica.tElapsed;' ])
    
        % class probabilities condensed to 4 categories
        ica.classProb4_classLabs       = v.ica.iclabel.classProb4.classLabs; % class labels for each of the four component categories
        eval([ 'ica.m_classProb4_mean(' pstr ',:)       = v.ica.iclabel.classProb4.mean;']) % mean probability that a component belongs to each category
        eval([ 'ica.m_classProb4_medi(' pstr ',:)       = v.ica.iclabel.classProb4.medi;']) % median probability that a component belongs to each category
        eval([ 'ica.m_classProb4_catNumProp(' pstr ',:) = v.ica.iclabel.classProb4.catNumProp;']) % percentage of components belonging to each category
        eval([ 'ica.m_classProb4_catVarProp(' pstr ',:) = v.ica.iclabel.classProb4.catVarProp;']) % percentage of total variance accounted for by each category
    
        % identifiability summary metrics
        eval([ 'ica.m_identiDiff_meanDiff('     pstr ')  = v.ica.iclabel.identiDiff.meanDiff;'])  % (diff prob between most likely category and unknown (other) categories)
        eval([ 'ica.m_identiDiff_propNumKnown(' pstr ')  = v.ica.iclabel.identiDiff.propNumKnown;'])  % number of components which are most likely to be a known category
        eval([ 'ica.m_identiDiff_propVarKnown(' pstr ')  = v.ica.iclabel.identiDiff.propVarKnown;'])  % percentage of total variance accounted for by known categories
        
        end
    end
        
        
    end % end file loop for 1st order metrics
    
    %% handle cleaning parameters which were varied
    % get a string which allows indexing each output matrix by the original data value
    if npars > 0
        pstr = repmat('1,',1,npars);                
        pstr(end) = []; % remove trailing comma     
    else
        pstr = '1';
    end
    
    %% compute second order metrics
    
    % fft
    if exist('fft','var')
        % compute proportions of bin amplitude relative to uncleaned dataset
        eval([ 'fft.m_binAmp_prop          = 100*(fft.m_binAmp  ./ fft.m_binAmp(' pstr ',:,:));']) 
        eval([ 'fft.m_binAmp_avgChan_prop  = 100*(fft.m_binAmp_avgChan ./ fft.m_binAmp_avgChan(' pstr ',:));']) 
    end
    
    % sw
    if exist('sw','var')
        eval([ 'sw.m_nwaves_rel             = sw.m_nwaves - sw.m_nwaves(' pstr ');']) % number of slow waves relative to the number found in the uncleaned data
        eval([ 'sw.m_amp_mean_prop   = 100*(sw.m_amp_mean / sw.m_amp_mean(' pstr '));'])
        eval([ 'sw.m_amp_median_prop = 100*(sw.m_amp_median / sw.m_amp_median(' pstr '));'])
        eval([ 'sw.m_tval_median_rel = sw.m_tval_median - sw.m_tval_median(' pstr ');'])
    end
    
    % ica
    if exist('ica','var')
        eval(['ica.m_numIC_rel      = ica.m_numIC - ica.m_numIC(' pstr ');'])
        eval(['ica.m_tElapsed_prop  = 100*(ica.m_tElapsed / ica.m_tElapsed(' pstr '));'])
        eval(['ica.m_tElapsed_rel   = ica.m_tElapsed - ica.m_tElapsed(' pstr ');'])
        
        eval(['ica.m_classProb4_catNumProp_rel = ica.m_classProb4_catNumProp - ica.m_classProb4_catNumProp(' pstr ',:);'])
        eval(['ica.m_classProb4_catVarProp_rel = ica.m_classProb4_catVarProp - ica.m_classProb4_catVarProp(' pstr ',:);'])
        
        eval(['ica.m_identiDiff_meanDiff_rel     = ica.m_identiDiff_meanDiff     - ica.m_identiDiff_meanDiff(' pstr ');'])
        eval(['ica.m_identiDiff_propNumKnown_rel = ica.m_identiDiff_propNumKnown - ica.m_identiDiff_propNumKnown(' pstr ');'])
        eval(['ica.m_identiDiff_propVarKnown_rel = ica.m_identiDiff_propVarKnown - ica.m_identiDiff_propVarKnown(' pstr ');'])
    end
    
    
    %% store outputs
    temp.pars = pars; 
    temp.asr = asr;
    if exist('fft',   'var'),     temp.fft    = fft;    end
    if exist('sw',    'var'),     temp.sw     = sw;     end
    if exist('ica',   'var'),     temp.ica    = ica;    end
    if nstages > 1,     temp.stage  = stage;  end
    
    %% if multiple parameter dimensions, recursively set all repeats of the orig data in output structure to be equal to the first element in the output matrix (original data)
    if npars > 1
        temp = rec_fixOrigVals(temp, npars);
    end

    %% store output for this stage
    valid_merged(st) = temp;


end % STAGE LOOP

%% save merged validation structure structure in each dataset
if ~isempty(cfg.loadPath)
    ctemp = []; 
    ctemp.headerOnly = true;
    ctemp.savePath   = cfg.loadPath;
    for f = 1 : length(cfg.loadNames)
        % load, add validation struct and resave
        EEG = pop_loadset('filepath', cfg.loadPath, 'filename', cfg.loadNames{f}, 'loadmode', 'info');
        EEG.etc.dusk2dawn.valid_merged = valid_merged;
        ctemp.saveName = cfg.loadNames{f};
        d2d_saveDataset(EEG, ctemp);
    end
else % if running fully in RAM
    EEG = evalin("caller",'EEG');
    EEG.etc.dusk2dawn.valid_merged = valid_merged;
    for f = 1 : length(cfg.loadNames)-1
        EEG.etc.dusk2dawn.EEG_cleaned.EEG{f}.etc.dusk2dawn.valid_merged = valid_merged;
    end
    assignin("caller","EEG",EEG);
end

fprintf('... finished merging validation metrics across ASR thresholds\n')

end % END FUNCTION

%% subfunction - recursively set the repeats of the metrics for original data to the same as the first element
function out = rec_fixOrigVals(in, npars)
   % get fields
   flds = fields(in);
   
   % loop through fields
   for f = 1:length(flds)
       
       % get field
       field = in.(flds{f});
       
       % if struct, apply function recursively
       if isstruct(field)
           field = rec_fixOrigVals(field, npars);
           out.(flds{f}) = field;
       else % otherwise check if field is numerical & a metric which needs to be duplicated, and if so then get all indices corresponding to original data 
%              and set these to be equal to the original data metric
            if isnumeric(field) && startsWith( flds{f} , 'm_' )
                
                   % construct indices referring to the element containing the original data value for this metric
                   extraDims = ndims(field)-npars;
                   indstr_rhs =  [ repmat({'1'},1,npars),  repmat({':'},1, extraDims) ];
                   indstr_rhs = strjoin(indstr_rhs,',');
                   for p = 1:npars
                       % construct indices spanning each corner formed by one parameter dimension and all the other parameters
                       indstr_lhs = repmat({':'},1,ndims(field));
                       indstr_lhs{p} = '1';
                       indstr_lhs = strjoin(indstr_lhs,','); 
                       % set those elements of field to be equal to the original data value for this metric
                       eval([ 'sizeNeeded = size(field('  indstr_lhs '));' ]) % get number of copies needed for assignment
                       if npars == 3 && length(sizeNeeded) == (npars-1) % if last dim is singleton, then size(field(indstr_lhs)), will not have the last "1", so add this back in
                            sizeNeeded = [sizeNeeded, 1]; %? this happens when there are 3 varied parameters and you reach p = 3
                       end
                       copystr  = num2str([ sizeNeeded(1:npars)  ones(1,extraDims) ], '%d,'); copystr(end) = []; % copystr is ones at the end of any extra dims (extra rel. to the varied parameter dimensions)
                       eval([ 'field(' indstr_lhs ') = repmat(  field(' indstr_rhs ') , ' copystr  ' ) ;'  ])             %  -> because you don't need to copy those dims which are already indexed as ":" by indstr_rhs
                   end % loop through parameters
                   
            end % if numeric
       end % if field is structure
       
       % set field
       out.(flds{f}) = field;
   end
end
