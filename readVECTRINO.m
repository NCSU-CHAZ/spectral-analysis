function [stVECT] = readVECTRINO(sVECTfldr)
% 
% Purpose: To read vectrino data (from multiple instruments) into a 
%          structure. Does not allow for split files. Requires time stamp
%          output option to be chosen within Nortek post-processing
%          software. Although there are four beams, we ignore the second z2
%          beam and associated QC criteria (for the side looker, the z2 
%          velocity is zero). As usual, Nortek sucks in that it does not
%          report the start or end time with millisecond precision in the 
%          header file; therefore, there will be a lag if compared to a 
%          non-synced instrument. However, the Vectrino profiler (.mat)
%          does report fractional seconds! We interpret the measurement 
%          time as the start date/first record + the instrument timestamp.
%          Lastly, the time of first/last measurement in the header file of
%          slaves (if synced with a master) will be incorrect: use the
%          master dt.
%
% Inputs:
%       - sVECTfldr: string of file name containing Vectrino data 
% 
% See also: readVECT.m
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       5/9/17          KA                  Original code 
%
%---------------------------------preamble--------------------------------%

disp('-----------------------------------------------------------')
disp('----------------------readVECTRINO-------------------------')
disp('-----------------------------------------------------------')               

% preallocate structure
stVECT = struct('velX',{}, 'velY',{}, 'velZ1',{}, ...
                'dt',{}, 'S2Nx',{}, 'S2Ny',{}, 'S2Nz1',{}, ...
                'corrX',{}, 'corrY',{}, 'corrZ1',{});

%---------------------------------read data-------------------------------%

% read filenames 2 structure        
stDatFiles = dir(fullfile(sVECTfldr, '*.dat')); % read filename 2 structure
stHdrFiles = dir(fullfile(sVECTfldr, '*.hdr')); 
stMatFiles = dir(fullfile(sVECTfldr, '*.mat')); 

for i=1 : length(stDatFiles)
    % read header data (start/end time and sampling rate)
    fid1       = fopen(fullfile(sVECTfldr, stHdrFiles(i).name));
    cDT        = textscan(fid1, '%*s %*s %*s %*s %s %s %s', 2, ...
                          'headerlines', 4); % read datetime data
    tempDT     = datetime(strcat(cDT{1}, {' '}, cDT{2}, {' '}, cDT{3}),...
                          'InputFormat', 'MM/dd/yyyy hh:mm:ss a');
    
    % read velocity and QC data 
    fid2      = fopen(fullfile(sVECTfldr, stDatFiles(i).name));
    cVECT     = textscan(fid2, ...   % read data
         '%f %*f %*f %f %f %f %*f %*f %*f %*f %*f %f %f %f %*f %f %f %f %*[^\n]'); 
    
    fclose('all');
        
    % save velocity data
    stVECT(i).velX  = [cVECT{2}]; % (m/s) 
    stVECT(i).velY  = [cVECT{3}]; % (m/s) 
    stVECT(i).velZ1 = [cVECT{4}]; % (m/s) 
    
    % save signal to noise ratio 
    stVECT(i).S2Nx  = [cVECT{5}]; % (dB)
    stVECT(i).S2Ny  = [cVECT{6}]; % (dB)
    stVECT(i).S2Nz1 = [cVECT{7}]; % (dB)
    
    % save correlation 
    stVECT(i).corrX = [cVECT{8}]; % (%)
    stVECT(i).corrY = [cVECT{9}]; % (%)
    stVECT(i).corrZ1= [cVECT{10}];% (%)   
    
    % create datetime array
    stVECT(i).dt    = tempDT(1) + seconds(cVECT{1});

end

for i=1 : length(stMatFiles)
    
    % load .mat file
    fid3 = fullfile(sVECTfldr, stMatFiles(i).name);
    load(fid3) 
    
    % save velocity data (profiles)
    stVECT(i).velX = Data.Profiles_VelX;  % (m/s) 
    stVECT(i).velY = Data.Profiles_VelY;  % (m/s) 
    stVECT(i).velZ1= Data.Profiles_VelZ1; % (m/s) 
    
    % save signal to noise ratio 
    stVECT(i).S2Nx      = Data.Profiles_SNRBeam1; % (dB)
    stVECT(i).S2Ny      = Data.Profiles_SNRBeam2; % (dB)
    stVECT(i).S2Nz1     = Data.Profiles_SNRBeam3; % (dB)
    
    % save correlation 
    stVECT(i).corrX     = Data.Profiles_CorBeam1; % (%)
    stVECT(i).corrY     = Data.Profiles_CorBeam2; % (%)
    stVECT(i).corrZ1    = Data.Profiles_CorBeam3; % (%)
    
    % create datetime array
%     tempDT = datetime(Data.Profiles_firstRecord, 'ConvertFrom', ...
%                       'posixtime', 'TimeZone', '-05:00');
    tempDT = datetime(Data.Profiles_firstRecord, 'ConvertFrom', ...
                      'posixtime') - hours(5);
    stVECT(i).dt = tempDT(1) + seconds(Data.Profiles_TimeStamp);
    
    close all
end

disp(' ')
disp('~~~~~~~~~~~~~~~~~~Finished reading Vectrino data~~~~~~~~~~~~~~~~~~~') 
disp(' ')

end