function [stVECT] = readVECT(sVECTfldr, bSEN)
% 
% Purpose: To read Vector data (from one instrument) into a structure.
%          Allows for split files.
%
% Note on time conversion: 
% I can only get the time stamp to match up correctly to the number of
% elements in the data arrays when 1) I convert the data from the post-
% processing software using an even time period where I know the Vector is
% in the water (it does not work on the full time series - this might 
% have to do with known errors when the Vector is in air), and 2) I add a 
% second to the end of the datetime Vector. See detailed code comments 
% below. Remember, you can use second() to check the precision in ms on the 
% datetime value! AND you can use xcorr to check the correlation and lag
% between signals.
%
% Inputs:
%       - sVECTfldr: string of file name containing Vector data ('Vector')
%       - bSEN: boolean for reading .SEN data
% 
% Record of revisions:
%       Date            Programmer    Description of Change
%       =========================================================
%       4/18/16         KA            Original code 
%       6/15/16         KA            Changed X,Y,Z to E,N,Up and modified
%                                     code to concatenate split files
%       9/15/17         KA            Added capabilities to read .SEN files
%      
%---------------------------------preamble--------------------------------%

disp('-----------------------------------------------------------')
disp('------------------------readVECT---------------------------')
disp('-----------------------------------------------------------')               

% preallocate structure
stTEMP = struct('velocityE',{}, 'velocityN',{}, 'velocityUp',{}, ...
                'dt',{}, 'pressure',{}, 'S2Neast',{}, 'S2Nnorth',{}, ...
                'S2Nup',{}, 'corrEast',{}, 'corrNorth',{}, 'corrUp',{}, ...
                'OBS1',{}, 'OBS2',{}, 'heading',{}, 'pitch',{}, 'roll',{}, ...
                'dtSEN',{});
stVECT = struct('velocityE',{}, 'velocityN',{}, 'velocityUp',{}, ...
                'dt',{}, 'pressure',{}, 'S2Neast',{}, 'S2Nnorth',{}, ...
                'S2Nup',{}, 'corrEast',{}, 'corrNorth',{}, 'corrUp',{}, ...
                'OBS1',{}, 'OBS2',{}, 'heading',{}, 'pitch',{}, 'roll',{}, ...
                'dtSEN',{});
            
%----------------------------read Vector data-----------------------------%

% load vector .DAT data into 2D array ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stDatFiles = dir(fullfile(sVECTfldr, '*.dat')); % read filename 2 structure
stHdrFiles = dir(fullfile(sVECTfldr, '*.hdr')); 
stSenFiles = dir(fullfile(sVECTfldr, '*.sen')); 

for i=1 : length(stDatFiles)
    % read .HDR data (start/end time and sampling rate)
    fid1       = fopen(fullfile(sVECTfldr, stHdrFiles(i).name));
    cDT        = textscan(fid1, '%*s %*s %*s %*s %s %s %s', 2, ...
                          'headerlines', 6); % read datetime data
    tempDT(:,i)= datetime(strcat(cDT{1}, {' '}, cDT{2}, {' '}, cDT{3}),...
                          'InputFormat', 'MM/dd/yyyy hh:mm:ss a');
    cHz        = textscan(fid1, '%*s %*s %f %*s', 1, ...
                        'headerlines', 4); % read sampling rate
    nHz        = [cHz{1}]; 
    
    % read .DAT data (x velocity and signal-2-noise)
    fid2      = fopen(fullfile(sVECTfldr, stDatFiles(i).name));
    cVect     = textscan(fid2, ...   % read data
        '%*f %*f %f %f %f %*f %*f %*f %f %f %f %f %f %f %f %f %f %*[^\n]'); 
    
    % save velocity data
    stTEMP(i).velocityE  = [cVect{1}]; % (m/s) 
    stTEMP(i).velocityN  = [cVect{2}]; % (m/s) 
    stTEMP(i).velocityUp = [cVect{3}]; % (m/s) 
    
    % save signal to noise ratio 
    stTEMP(i).S2Neast    = [cVect{4}]; % (dB)
    stTEMP(i).S2Nnorth   = [cVect{5}]; % (dB)
    stTEMP(i).S2Nup      = [cVect{6}]; % (dB)
    
    % save correlation 
    stTEMP(i).corrEast   = [cVect{7}]; % (dB)
    stTEMP(i).corrNorth  = [cVect{8}]; % (dB)
    stTEMP(i).corrUp     = [cVect{9}]; % (dB)   
    
    % save pressure
    stTEMP(i).pressure   = [cVect{10}]; % (dbar)
    
    % save OBS data
    stTEMP(i).OBS1       = [cVect{11}]; % (NTU)
    stTEMP(i).OBS2       = [cVect{12}]; % (NTU)
    
    if bSEN
        % read .SEN data (heading, pitch, roll)
        fid3      = fopen(fullfile(sVECTfldr, stSenFiles(i).name));
        cSen      = textscan(fid3, ...   % read data
            '%f %f %f %f %f %f %*f %*f %*f %*f %f %f %f %*[^\n]');
    
        % save heading, pitch, roll data and datetime collected
        stTEMP(i).dtSEN   = datetime(cSen{3}, cSen{1}, cSen{2}, cSen{4}, ...
                                     cSen{5}, cSen{6});
        stTEMP(i).heading = [cSen{7}]; % (degrees) 
        stTEMP(i).pitch   = [cSen{8}]; % (degrees)
        stTEMP(i).roll    = [cSen{9}]; % (degrees) 
    end
    
    fclose('all');
end

% concatenate each field into a single matrix
stVECT(1).pressure    = cat(1,stTEMP.pressure);
stVECT(1).OBS1        = cat(1,stTEMP.OBS1);     
stVECT(1).OBS2        = cat(1,stTEMP.OBS2);
stVECT(1).velocityE   = cat(1,stTEMP.velocityE);
stVECT(1).velocityN   = cat(1,stTEMP.velocityN);
stVECT(1).velocityUp  = cat(1,stTEMP.velocityUp);
stVECT(1).S2Neast     = cat(1,stTEMP.S2Neast);
stVECT(1).S2Nnorth    = cat(1,stTEMP.S2Nnorth);
stVECT(1).S2Nup       = cat(1,stTEMP.S2Nup);
stVECT(1).corrEast    = cat(1,stTEMP.corrEast);
stVECT(1).corrNorth   = cat(1,stTEMP.corrNorth);
stVECT(1).corrUp      = cat(1,stTEMP.corrUp);
stVECT(1).heading     = cat(1,stTEMP.heading);
stVECT(1).pitch       = cat(1,stTEMP.pitch);
stVECT(1).roll        = cat(1,stTEMP.roll);
if bSEN
    stVECT(1).dtSEN       = cat(1,stTEMP.dtSEN);
end
%-------------------------create datetime array---------------------------%

% If you tell the Vector to start at 12:00:00, the first velocity sample 
% begins at 12:00:02.0 (see Research/BARis/ADV_ADCP/VectorTimeConversion). 
% Therefore, if you are sampling at 16 Hz, the first velocity sample is 
% complete at 12:00:02.0625. Since the Vector is pinging continuously 
% during each sample, the best time for the sample would be the midpoint of 
% the measurement. So the time for the first 16 Hz sample would be 
% 12:00.02.0312: 
%     stVECT(1).dt = (tempDT(1,1) + seconds(2 + 1/32) : seconds(1/16) : ...
%                     tempDT(2,length(stDatFiles)) + seconds(3))';
% BUT, this is super annoying for comparing time series so I decided to
% assign the timestamp to the end of the measurement period:
stVECT(1).dt = (tempDT(1,1) + seconds(2 + 1/nHz) : seconds(1/nHz) : ...
                tempDT(2,length(stDatFiles)))';
            
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~Finished reading Vector data~~~~~~~~~~~~~~~~~~~~') 
disp(' ')

return