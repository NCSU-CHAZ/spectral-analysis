function [stATM] = readATM(sATMfldr)
% 
% Purpose: To read atmospheric pressure data from a HOBO or from a NOAA 
% meteorological station. NOAA files must be modified such that the 2nd 
% column is datetime and the 3rd column is barometric pressure; add one row
% to the top; convert mb to kPa (divide by ten). Do not store any other 
% files in the same directory. For MacOs 
% users, delete .DS_STORE in terminal.
%
% Inputs:
%       - sATMfldr: string of filename containing barometric pressure data
% 
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       9/5/17          KAnarde             Original code      
%
%% ---------------------------------preamble---------------------------- %%

disp('-----------------------------------------------------------')
disp('-------------------------readATM---------------------------')
disp('-----------------------------------------------------------')               

% preallocate structure which is organized by POD # (again, needs to be
% part of the data filename)

stATM = struct('baropress',{}, 'dt',{});

%% ----------------------read barometric pressure data------------------ %%
% read filenames 2 structure
stATMfiles = dir(sATMfldr);
stATMfiles = stATMfiles(~[stATMfiles.isdir]); % keep files, not directories
                    
for i = 1 : numel(stATMfiles)            

    fid = fopen(fullfile(sATMfldr, stATMfiles(i).name));
    fgetl(fid); fgetl(fid);               % ignore first 2 lines, junk
    cATM = textscan(fid, '%*s %s %f %*f', 'delimiter', ',');
    fclose(fid);

    % save date and time to datetime variable
    dtATMp           = cATM{1};
    stATM(i).dt      = datetime(dtATMp,'InputFormat', ...
                                'MM/dd/yy hh:mm:ss a');

    % save pressure (barometric pressure [kPa])
    stATM(i).baropress= cATM{2};
end        
        
disp(' ')
disp('~~~~~~~~~~~~~~~~~Finished reading ATM data~~~~~~~~~~~~~~~~~~~~')
disp(' ')
return