function [stPOD] = readPOD(sPSfldr, sTCMfldr)
% 
% Purpose: To read pressure sensor (HOBO and RBR) and TCM data into a 
% structure that is organized by POD # (integers). Requires that the data 
% be arranged by folder and that the files within the folder are ordered 
% correctly by increasing POD #. Do not store any other files in 
% the same directory. For MacOs users, delete .DS_STORE in terminal.
%
% Inputs:
%       - sPSfldr: string of filename containing pressure sensor data('PS')
%       - sTCMfldr: string of filename containing TCM data ('TCM')
% 
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       4/18/16         KAnarde             Original code 
%       6/13/16         KA                  Modified for HOBOs 
%       2/20/17         KA                  Modified for RSKtools       
%
%% ---------------------------------preamble---------------------------- %%

disp('-----------------------------------------------------------')
disp('-------------------------readPOD---------------------------')
disp('-----------------------------------------------------------')               

% preallocate structure which is organized by POD # (again, needs to be
% part of the data filename)

stPOD = struct('pressure',{}, 'seapress',{}, 'depth',{},'baropress', {},...
    'PSdt',{}, 'PSunits',{}, 'PSid',{}, 'PStemp',{}, 'speed',{}, ...
    'bearing',{}, 'roll', {}, 'pitch', {}, 'yaw', {}, 'velocityN',{}, ...
    'velocityE',{}, 'TCMdt',{}, 'TCMunits',{}, 'Ax', {}, 'Ay', {}, ...
    'Az', {}, 'Mx', {}, 'My', {}, 'Mz', {});

%% ------------------------read pressure sensor data-------------------- %%
% read filenames 2 structure
stPSfiles = dir(sPSfldr);
stPSfiles = stPSfiles(~[stPSfiles.isdir]); % keep files, not directories
                    
for i = 1 : numel(stPSfiles)            
    [~,~,sExt] = fileparts(stPSfiles(i).name);  % get extension
    
    % HOBOs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if strcmp(sExt,'.csv')
        fid = fopen(fullfile(sPSfldr, stPSfiles(i).name));
        fgetl(fid);                      % ignore first line, junk
        stPOD(i).PSunits = fgetl(fid);   % read 2nd header
        cHOBO = textscan(fid, '%*s %s %f %f %f %f','delimiter', ',');
        fclose(fid);
        
        % save date and time to datetime variable
        dtHOBOp           = cHOBO{1};
        stPOD(i).PSdt     = datetime(dtHOBOp,'InputFormat', ...
                                     'MM/dd/yy hh:mm:ss a ');
        
        % save pressure (absolute [kPa], temperature [C], barometric
        % pressure [kPa], depth [m])
        stPOD(i).pressure = cHOBO{2};
        stPOD(i).PStemp   = cHOBO{3};
        stPOD(i).baropress= cHOBO{4};
        stPOD(i).depth    = cHOBO{5};
        
        % REMOVE NaNs!
        % first, NaNs exist because of the offset with AtmPressure readings 
        % so we must interpolate the AtmPressure 
        iNaN = isnan(stPOD(i).baropress);
        cnt  = 1:numel(stPOD(i).baropress);
        stPOD(i).baropress(iNaN) = interp1(cnt(~iNaN), ...
                                       stPOD(i).baropress(~iNaN), ...
                                       cnt(iNaN));
        % next, remove all elements associated with the ~NaNs in pressure
        % and depth columns, i.e. we only want to keep the new AtmPressure
        % elements which correspond with the indices of the other data sets
        iNaN = isnan(stPOD(i).pressure);
        stPOD(i).pressure = stPOD(i).pressure(~iNaN);
        stPOD(i).PStemp   = stPOD(i).PStemp(~iNaN);
        stPOD(i).baropress= stPOD(i).baropress(~iNaN);
        stPOD(i).depth    = stPOD(i).depth(~iNaN);  
        stPOD(i).PSdt     = stPOD(i).PSdt(~iNaN); 
    
    % RBRs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(sExt,'.rsk') 
        
        % open .rsk file using RSKtools (downloaded from RBR); this 
        % generates a structure with the variable names and a thumbnail of 
        % the data
        RBR = RSKopen(fullfile(sPSfldr,stPSfiles(i).name));
        
        % read all of the data from the database to the disk (RSKreaddata 
        % can also just read portions of the data w/ a start and end time);
        % for large files this can take upwards of 10 minutes
        RBR = RSKreaddata(RBR);
        
        % save pressure (dbar), datetime, and units
        stPOD(i).pressure = RBR.data.values(:,1);
        stPOD(i).PSdt     = datetime(RBR.data.tstamp,'ConvertFrom', ...
                                     'datenum');
        stPOD(i).PSunits  = RBR.channels.units;
        stPOD(i).PSid     = RBR.instruments.serialID;
        
        % clear structure 
        clear RBR
        close all
    end
end

%% -------------------------------read TCM data------------------------- %%

% first process the current data (_CR.txt) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stFiles = dir(fullfile(sTCMfldr, '*CR.txt')); % read filenames 2 structure

for i=1 : numel(stFiles)
    sFileName = stFiles(i).name;    % retrieve filename from structure
    fid  = fopen(fullfile(sTCMfldr, sFileName));
    stPOD(i).TCMunits = fgetl(fid); % read 1st line of file
    % use textscan because of weird datetime text
    cTCM = textscan(fid, '%s %f %f %f %f','delimiter', ',');
    fclose(fid);

    % save datetime, speed (cm/s), bearing (degrees), vel-N (cm/s), and
    % vel-E (cm/s) 
    stPOD(i).TCMdt     = datetime([cTCM{1}],'InputFormat', ... 
                                  'uuuu-MM-dd''T''HH:mm:ss.SSS'); 
    stPOD(i).speed     = cTCM{2};   
    stPOD(i).bearing   = cTCM{3};
    stPOD(i).velocityN = cTCM{4};
    stPOD(i).velocityE = cTCM{5};
end

% next process the roll/pitch/yaw data (_RPY.txt) ~~~~~~~~~~~~~~~~~~~~~~~~~
stFiles = dir(fullfile(sTCMfldr, '*RPY.txt')); % read filenames 2 structure

for i=1 : numel(stFiles)
    sFileName = stFiles(i).name;    % retrieve filename from structure
    fid  = fopen(fullfile(sTCMfldr, sFileName));
    fgetl(fid);                     % skip units
    % use textscan because of weird date/time text
    cTCMrpy = textscan(fid, '%*s %f %f %f %*f %*f %*f','delimiter', ',');
    fclose(fid);

    % save RPY data (degrees)
    stPOD(i).roll  = cTCMrpy{1};   
    stPOD(i).pitch = cTCMrpy{2};
    stPOD(i).yaw   = cTCMrpy{3};
end

% last, process the accelerometer and magnetometer data (_MA.txt)
stFiles = dir(fullfile(sTCMfldr, '*MA.txt')); % read filenames 2 structure

for i=1 : numel(stFiles)
    sFileName = stFiles(i).name;    % retrieve filename from structure
    fid  = fopen(fullfile(sTCMfldr, sFileName));
    fgetl(fid); % skip units, we know = [g] 
    cTCMma = textscan(fid, '%s %f %f %f %f %f %f','delimiter', ',');
    fclose(fid);

    % save datetime; accelerometer data Ax, Ay, Az (g); and magnetometer
    % data Mx, My, Mx (mg)
    stPOD(i).TCMdtA = datetime([cTCMma{1}],'InputFormat', ...
                                           'uuuu-MM-dd''T''HH:mm:ss.SSS'); 
    stPOD(i).Ax = cTCMma{2};   
    stPOD(i).Ay = cTCMma{3};
    stPOD(i).Az = cTCMma{4};
    stPOD(i).Mx = cTCMma{5};   
    stPOD(i).My = cTCMma{6};
    stPOD(i).Mz = cTCMma{7};
end

disp(' ')
disp('~~~~~~~~~~~~~~~~~Finished reading POD data~~~~~~~~~~~~~~~~~~~~')
disp(' ')
return