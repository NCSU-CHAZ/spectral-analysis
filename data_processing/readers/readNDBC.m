function [stNDBC] = readNDBC(sNDBCfldr)
% 
% Purpose: To read NDBC data from gauge 42019 (82 m water depth). Do not 
% store other files in the same directory. For MacOs users, delete 
% .DS_STORE in terminal. Note that the reading is hard-coded such that
% alpha1 is read first, followed by alpha2, r1, r2, and psd.
%
% Inputs:
%       - sNDBCfldr: string of filename containing NDBC spectral wave data
% 
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       01/5/18         KAnarde             Original code    
%% ---------------------------------preamble---------------------------- %%

disp('-----------------------------------------------------------')
disp('-------------------------readNDBC--------------------------')
disp('-----------------------------------------------------------')               

% preallocate structure 
stNDBC = struct('dt', {}, 'a1',{}, 'a2', {}, 'r1', {}, 'r2', {}, ...
                'psd', {}, 'f', {});

% brute force these things because I'm too lazy to code something
nNumberCols = 52;

%% -------------------------------read data----------------------------- %%
% read filenames 2 structure
stNDBCfiles = dir(sNDBCfldr);
stNDBCfiles = stNDBCfiles(~[stNDBCfiles.isdir]); % don't keep directories
                    
for i = 1 : numel(stNDBCfiles)   
           
    [~,~,sExt] = fileparts(stNDBCfiles(i).name);  % get extension
    
    % read text files only
    if strcmp(sExt,'.txt')
        fid = fopen(fullfile(sNDBCfldr, stNDBCfiles(i).name));
        % header = frequency bins
        format = ['%*s %*s %*s %*s %*s' repmat('%f', [1 nNumberCols-5])];
        stNDBC{1,1}.f = cell2mat(textscan(fid, format, 1))';
        % first five columns are YY MM DD HH MM
        format = repmat('%f', [1 nNumberCols]);
        cDATA = textscan(fid, format);
        fclose(fid);
        
        % save date and time to datetime variable 
        stNDBC{1,1}.dt = datetime([cDATA{1}, cDATA{2}, cDATA{3}, cDATA{4}, ...
                              cDATA{5}, zeros(size(cDATA{1}))])';
        
        % save the rest of the data  
        if i == 1
            stNDBC{1,1}.a1  = cell2mat(cDATA(6:end))';
        elseif i == 2
            stNDBC{1,1}.a2  = cell2mat(cDATA(6:end))';
        elseif i == 3
            stNDBC{1,1}.r1  = cell2mat(cDATA(6:end))';
        elseif i == 4
            stNDBC{1,1}.r2  = cell2mat(cDATA(6:end))';
        elseif i == 5
            stNDBC{1,1}.psd = cell2mat(cDATA(6:end))';
        end
    end
end        
        
disp(' ')
disp('~~~~~~~~~~~~~~~~~Finished reading NDBC data~~~~~~~~~~~~~~~~~~~~')
disp(' ')
return