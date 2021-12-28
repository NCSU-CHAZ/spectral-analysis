% ------------------------runCOASTRR_DataReport---------------------------%
% Purpose: This run script creates netcdf files for reporting on
% DesignSafe; note that the file structure is based on ADCIRC output files
%
% Record of revisions:
%       Date            Programmer          Description of Change
%       =========================================================
%       09/03/19        KAnarde             Original code
%
%% PREAMBLE

% change to the data directory on external hard-drive 
cd /Volumes/Katherine_4TB/BARis/COASTRR/HARVEY/DATA/POSTSTORM/
load('RRU_Reporting')
load('ADV_Reporting')
load('VECT')

%% create netcdf file for PT-1

% create variables
iStart = 2764801; iEnd = 9814282; % this is the windowed data set
nccreate('PT1.nc','pressure', 'Dimensions', {'time', length(stRRU.PT(1).seapress(iStart:iEnd))})
nccreate('PT1.nc','dt', 'Dimensions', {'time', length(stRRU.PT(1).dt(iStart:iEnd))})

% write data to file
ncwrite('PT1.nc', 'pressure', stRRU.PT(1).seapress(iStart:iEnd))
ncwrite('PT1.nc', 'dt', datenum(stRRU.PT(1).dt(iStart:iEnd)+hours(5))) % convert to UTC from LDT

% create local attributes (units and description)
ncwriteatt('PT1.nc','dt', 'units', 'datenum');
ncwriteatt('PT1.nc','dt', 'description', 'time relative to UTC (LDT + 5 hours)');
ncwriteatt('PT1.nc','pressure', 'units', 'dbar');
ncwriteatt('PT1.nc','pressure', 'description', 'gauge pressure, meaning pressure corrected for atmospheric/barometric pressure; barometric pressure data was downloaded from https://tidesandcurrents.noaa.gov for NOAA Station 8773146');

% testing
tmpDT = ncread('PT1.nc', 'dt');
tmpP = ncread('PT1.nc', 'pressure');
figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP)
% figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP, ...
%              stRRU.PT(1).dt, stRRU.PT(1).seapress)

% create global attributes (sensor description)
ncwriteatt('PT1.nc','/','description', 'Time series of gauge pressure recorded by a pressure transducer housed within a shallow PVC well located in the backshore at Matagorda Peninsula during Hurricane Harvey.');
ncwriteatt('PT2.nc','/','comments', 'See PT-1.pdf for depiction of instrument deployment setup and position relative to storm deposits.');
ncwriteatt('PT1.nc','/','creation_date', datestr(now));
ncwriteatt('PT1.nc','/','modification_date', 'none');
ncwriteatt('PT1.nc','/','institution', 'Texas A&M University, Department of Ocean Engineering');
ncwriteatt('PT1.nc','/','host', 'DesignSafe-CI');
ncwriteatt('PT1.nc','/','contact', 'kanarde@tamu.edu or figlusj@tamu.edu');
ncwriteatt('PT1.nc','/','sensor', 'RBR Dwave 41410');
ncwriteatt('PT1.nc','/','PT_samp_frequency', '16 Hz');
ncwriteatt('PT1.nc','/','PT_latitude', '28.612403097024846');
ncwriteatt('PT1.nc','/','PT_longitude', '-95.93936769918335');
ncwriteatt('PT1.nc','/','PT_elevation_prestorm', '0.8575 m NAVD88');
ncwriteatt('PT1.nc','/','PT_elevation_poststorm', '0.9315 m NAVD88');
ncwriteatt('PT1.nc','/','bed_elevation_prestorm', '0.8575 m NAVD88');
ncwriteatt('PT1.nc','/','bed_elevation_poststorm', '1.400 m NAVD88');
ncdisp('PT1.nc')

%% create netcdf file for PT-2

% create variables
iStart = 345601; iEnd = 1388881; % this is the windowed data set
nccreate('PT2.nc','pressure', 'Dimensions', {'time', length(stRRU.PT(2).seapress(iStart:iEnd))})
nccreate('PT2.nc','dt', 'Dimensions', {'time', length(stRRU.PT(2).dt(iStart:iEnd))})

% write data to file
ncwrite('PT2.nc', 'pressure', stRRU.PT(2).seapress(iStart:iEnd))
ncwrite('PT2.nc', 'dt', datenum(stRRU.PT(2).dt(iStart:iEnd)+hours(5))) % convert to UTC from LDT

% create local attributes (units and description)
ncwriteatt('PT2.nc','dt', 'units', 'datenum');
ncwriteatt('PT2.nc','dt', 'description', 'time relative to UTC (LDT + 5 hours)');
ncwriteatt('PT2.nc','pressure', 'units', 'dbar');
ncwriteatt('PT2.nc','pressure', 'description', 'gauge pressure, meaning pressure corrected for atmospheric/barometric pressure; barometric pressure data was downloaded from https://tidesandcurrents.noaa.gov for NOAA Station 8773146');

% testing
tmpDT = ncread('PT2.nc', 'dt');
tmpP = ncread('PT2.nc', 'pressure');
figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP)
%figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP, ...
%             stRRU.PT(2).dt, stRRU.PT(2).seapress)

% create global attributes (sensor description)
ncwriteatt('PT2.nc','/','description', 'Time series of gauge pressure recorded by a pressure transducer housed within a shallow PVC well located in the back-barrier at Matagorda Peninsula during Hurricane Harvey.');
ncwriteatt('PT2.nc','/','comments', 'See PT-2.pdf for depiction of instrument deployment setup and position relative to storm deposits.');
ncwriteatt('PT2.nc','/','creation_date', datestr(now));
ncwriteatt('PT2.nc','/','modification_date', 'none');
ncwriteatt('PT2.nc','/','institution', 'Texas A&M University, Department of Ocean Engineering');
ncwriteatt('PT2.nc','/','host', 'DesignSafe-CI');
ncwriteatt('PT2.nc','/','contact', 'kanarde@tamu.edu or figlusj@tamu.edu');
ncwriteatt('PT2.nc','/','sensor', 'RBR solo 77793');
ncwriteatt('PT2.nc','/','PT_samp_frequency', '2 Hz');
ncwriteatt('PT2.nc','/','PT_latitude', '28.616407057141966');
ncwriteatt('PT2.nc','/','PT_longitude', '-95.94081907817059');
ncwriteatt('PT2.nc','/','PT_elevation_prestorm', '0.222 m NAVD88');
ncwriteatt('PT2.nc','/','PT_elevation_poststorm', '0.260 m NAVD88');
ncwriteatt('PT2.nc','/','bed_elevation_prestorm', '0.222 m NAVD88');
ncwriteatt('PT2.nc','/','bed_elevation_poststorm', '0.384 m NAVD88');
ncdisp('PT2.nc')

%% create netcdf file for PT-3

% create variables
iStart = 345601; iEnd = 1301221; % this is the windowed data set
nccreate('PT3.nc','pressure', 'Dimensions', {'time', length(stRRU.PT(3).seapress(iStart:iEnd))})
nccreate('PT3.nc','dt', 'Dimensions', {'time', length(stRRU.PT(3).dt(iStart:iEnd))})

% write data to file
ncwrite('PT3.nc', 'pressure', stRRU.PT(3).seapress(iStart:iEnd))
ncwrite('PT3.nc', 'dt', datenum(stRRU.PT(3).dt(iStart:iEnd)+hours(5))) % convert to UTC from LDT

% create local attributes (units and description)
ncwriteatt('PT3.nc','dt', 'units', 'datenum');
ncwriteatt('PT3.nc','dt', 'description', 'time relative to UTC (LDT + 5 hours)');
ncwriteatt('PT3.nc','pressure', 'units', 'dbar');
ncwriteatt('PT3.nc','pressure', 'description', 'gauge pressure, meaning pressure corrected for atmospheric/barometric pressure; barometric pressure data was measured by a subaerial mounted PT at Follets Island (see ATM.nc).');

% testing
tmpDT = ncread('PT3.nc', 'dt');
tmpP = ncread('PT3.nc', 'pressure');
%figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP)
figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP, ...
             stRRU.PT(3).dt, stRRU.PT(3).seapress)

% create global attributes (sensor description)
ncwriteatt('PT3.nc','/','description', 'Time series of gauge pressure recorded by a pressure transducer housed within a shallow PVC well located in the back-barrier at Follets Island during Hurricane Harvey.');
ncwriteatt('PT3.nc','/','comments', 'See PT-3.pdf for depiction of instrument deployment setup.');
ncwriteatt('PT3.nc','/','creation_date', datestr(now));
ncwriteatt('PT3.nc','/','modification_date', 'none');
ncwriteatt('PT3.nc','/','institution', 'Texas A&M University, Department of Ocean Engineering');
ncwriteatt('PT3.nc','/','host', 'DesignSafe-CI');
ncwriteatt('PT3.nc','/','contact', 'kanarde@tamu.edu or figlusj@tamu.edu');
ncwriteatt('PT3.nc','/','sensor', 'RBR solo 77792');
ncwriteatt('PT3.nc','/','PT_samp_frequency', '2 Hz');
ncwriteatt('PT3.nc','/','PT_latitude', '29.03166495838045');
ncwriteatt('PT3.nc','/','PT_longitude', '-95.19516703557478');
ncwriteatt('PT3.nc','/','PT_elevation_prestorm', '0.219 m NAVD88');
ncwriteatt('PT3.nc','/','PT_elevation_poststorm', '0.190 m NAVD88');
ncwriteatt('PT3.nc','/','bed_elevation_prestorm', '0.219 m NAVD88');
ncwriteatt('PT3.nc','/','bed_elevation_poststorm', '0.171 m NAVD88');
ncdisp('PT3.nc')

%% create netcdf file for ATM

% create variables
iStart = 1; iEnd = length(stRRU.ATM(4).baropress); 
nccreate('ATM.nc','pressure', 'Dimensions', {'time', length(stRRU.ATM(4).baropress(iStart:iEnd))})
nccreate('ATM.nc','dt', 'Dimensions', {'time', length(stRRU.ATM(4).dt(iStart:iEnd))})

% write data to file
ncwrite('ATM.nc', 'pressure', stRRU.ATM(4).baropress(iStart:iEnd))
ncwrite('ATM.nc', 'dt', datenum(stRRU.ATM(4).dt(iStart:iEnd)+hours(5))) % convert to UTC from LDT

% create local attributes (units and description)
ncwriteatt('ATM.nc','dt', 'units', 'datenum');
ncwriteatt('ATM.nc','dt', 'description', 'time relative to UTC (LDT + 5 hours)');
ncwriteatt('ATM.nc','pressure', 'units', 'dbar');
ncwriteatt('ATM.nc','pressure', 'description', 'atmospheric/barometric pressure');

% testing
tmpDT = ncread('ATM.nc', 'dt');
tmpP = ncread('ATM.nc', 'pressure');
%figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP)
figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP, ...
             stRRU.ATM(4).dt, stRRU.ATM(4).baropress)

% create global attributes (sensor description)
ncwriteatt('ATM.nc','/','description', 'Time series of atmospheric pressure recorded by a subaerially mounted pressure transducer (on a telephone pole 20 feet above ground surface) at Follets Island during Hurricane Harvey.');
ncwriteatt('ATM.nc','/','comments', 'See ATM.pdf for depiction of instrument deployment setup.');
ncwriteatt('ATM.nc','/','creation_date', datestr(now));
ncwriteatt('ATM.nc','/','modification_date', 'none');
ncwriteatt('ATM.nc','/','institution', 'Texas A&M University, Department of Ocean Engineering');
ncwriteatt('ATM.nc','/','host', 'DesignSafe-CI');
ncwriteatt('ATM.nc','/','contact', 'kanarde@tamu.edu or figlusj@tamu.edu');
ncwriteatt('ATM.nc','/','sensor', 'Onset HOBO U20');
ncwriteatt('ATM.nc','/','PT_samp_period', '30 sec (0.033 Hz)');
ncwriteatt('ATM.nc','/','PT_latitude', '29.024395');
ncwriteatt('ATM.nc','/','PT_longitude', '-95.192153');
ncdisp('ATM.nc')

%% create netcdf file for ADV

% create variables
iStart = 345568; iEnd = length(stADV.DATA.seapress); % this is the windowed data set to account for barometric pressure
nccreate('ADV.nc','pressure', 'Dimensions', {'time', length(stADV.DATA.seapress(iStart:iEnd))})
nccreate('ADV.nc','dt', 'Dimensions', {'time', length(stADV.DATA.dt(iStart:iEnd))})
nccreate('ADV.nc','velocityE', 'Dimensions', {'time', length(stADV.DATA.velE(iStart:iEnd))})
nccreate('ADV.nc','velocityN', 'Dimensions', {'time', length(stADV.DATA.velN(iStart:iEnd))})
nccreate('ADV.nc','velocityU', 'Dimensions', {'time', length(stADV.DATA.velU(iStart:iEnd))})
nccreate('ADV.nc','snrE', 'Dimensions', {'time', length(stADV.QC.snrE(iStart:iEnd))})
nccreate('ADV.nc','snrN', 'Dimensions', {'time', length(stADV.QC.snrN(iStart:iEnd))})
nccreate('ADV.nc','snrU', 'Dimensions', {'time', length(stADV.QC.snrU(iStart:iEnd))})
nccreate('ADV.nc','corrE', 'Dimensions', {'time', length(stADV.QC.corrE(iStart:iEnd))})
nccreate('ADV.nc','corrN', 'Dimensions', {'time', length(stADV.QC.corrN(iStart:iEnd))})
nccreate('ADV.nc','corrU', 'Dimensions', {'time', length(stADV.QC.corrU(iStart:iEnd))})
nccreate('ADV.nc','heading', 'Dimensions', {'sensor_time', length(stADV.QC.heading)})
nccreate('ADV.nc','pitch', 'Dimensions', {'sensor_time', length(stADV.QC.pitch)})
nccreate('ADV.nc','roll', 'Dimensions', {'sensor_time', length(stADV.QC.roll)})
nccreate('ADV.nc','dtHPR', 'Dimensions', {'sensor_time', length(stVECT.dtSEN)})
nccreate('ADV.nc','ssc63cm', 'Dimensions', {'time', length(stVECT.BT2.OBS1c(iStart:iEnd))})
nccreate('ADV.nc','ssc35cm', 'Dimensions', {'time', length(stVECT.BT2.OBS2c(iStart:iEnd))})

% write data to file
ncwrite('ADV.nc', 'pressure', stADV.DATA.seapress(iStart:iEnd))
ncwrite('ADV.nc', 'dt', datenum(stADV.DATA.dt(iStart:iEnd)+hours(5))) % convert to UTC from LDT
ncwrite('ADV.nc', 'velocityE', stADV.DATA.velE(iStart:iEnd))
ncwrite('ADV.nc', 'velocityN', stADV.DATA.velN(iStart:iEnd))
ncwrite('ADV.nc', 'velocityU', stADV.DATA.velU(iStart:iEnd))
ncwrite('ADV.nc', 'snrE', stADV.QC.snrE(iStart:iEnd))
ncwrite('ADV.nc', 'snrN', stADV.QC.snrN(iStart:iEnd))
ncwrite('ADV.nc', 'snrU', stADV.QC.snrU(iStart:iEnd))
ncwrite('ADV.nc', 'corrE', stADV.QC.corrE(iStart:iEnd))
ncwrite('ADV.nc', 'corrN', stADV.QC.corrN(iStart:iEnd))
ncwrite('ADV.nc', 'corrU', stADV.QC.corrU(iStart:iEnd))
ncwrite('ADV.nc', 'heading', stADV.QC.heading)
ncwrite('ADV.nc', 'pitch', stADV.QC.pitch)
ncwrite('ADV.nc', 'roll', stADV.QC.roll)
ncwrite('ADV.nc', 'dtHPR', datenum(stVECT.dtSEN+hours(5))) % convert to UTC from LDT
ncwrite('ADV.nc', 'ssc63cm', stVECT.BT2.OBS1c(iStart:iEnd))
ncwrite('ADV.nc', 'ssc35cm', stVECT.BT2.OBS2c(iStart:iEnd))

% create local attributes (units and description)
ncwriteatt('ADV.nc','dt', 'units', 'datenum');
ncwriteatt('ADV.nc','dt', 'description', 'time relative to UTC (LDT + 5 hours)');
ncwriteatt('ADV.nc','pressure', 'units', 'dbar');
ncwriteatt('ADV.nc','pressure', 'description', 'gauge pressure, meaning pressure corrected for atmospheric/barometric pressure; barometric pressure data was measured by a subaerial mounted PT at Follets Island (see ATM.nc).');
ncwriteatt('ADV.nc','velocityE', 'units', 'm/s');
ncwriteatt('ADV.nc','velocityE', 'description', 'easting velocity');
ncwriteatt('ADV.nc','velocityN', 'units', 'm/s');
ncwriteatt('ADV.nc','velocityN', 'description', 'northing velocity');
ncwriteatt('ADV.nc','velocityU', 'units', 'm/s');
ncwriteatt('ADV.nc','velocityU', 'description', 'upward velocity');
ncwriteatt('ADV.nc','snrE', 'units', 'dB');
ncwriteatt('ADV.nc','snrE', 'description', 'easting signal-to-noise ratio');
ncwriteatt('ADV.nc','snrN', 'units', 'dB');
ncwriteatt('ADV.nc','snrN', 'description', 'northing signal-to-noise ratio');
ncwriteatt('ADV.nc','snrU', 'units', 'dB');
ncwriteatt('ADV.nc','snrU', 'description', 'upward signal-to-noise ratio');
ncwriteatt('ADV.nc','corrE', 'units', '%');
ncwriteatt('ADV.nc','corrE', 'description', 'easting percent correlation');
ncwriteatt('ADV.nc','corrN', 'units', '%');
ncwriteatt('ADV.nc','corrN', 'description', 'northing percent correlation');
ncwriteatt('ADV.nc','corrU', 'units', '%');
ncwriteatt('ADV.nc','corrU', 'description', 'upward percent correlation');
ncwriteatt('ADV.nc','heading', 'units', 'deg');
ncwriteatt('ADV.nc','heading', 'description', 'heading, note different time vector (dtHPR)');
ncwriteatt('ADV.nc','pitch', 'units', 'deg');
ncwriteatt('ADV.nc','pitch', 'description', 'pitch, note different time vector (dtHPR)');
ncwriteatt('ADV.nc','roll', 'units', 'deg');
ncwriteatt('ADV.nc','roll', 'description', 'roll, note different time vector (dtHPR)');
ncwriteatt('ADV.nc','dtHPR', 'units', 'datenum');
ncwriteatt('ADV.nc','dtHPR', 'description', 'sensor time (for heading, pitch, roll) relative to UTC (LDT + 5 hours)');
ncwriteatt('ADV.nc','ssc63cm', 'units', 'g/L');
ncwriteatt('ADV.nc','ssc63cm', 'description', 'suspended sediment concentration at 63 cm above the pre-storm bed elevation');
ncwriteatt('ADV.nc','ssc35cm', 'units', 'g/L');
ncwriteatt('ADV.nc','ssc35cm', 'description', 'suspended sediment concentration at 35 cm above the pre-storm bed elevation (concurrent with the measurement volume of the ADV)');

% testing
tmpDT = ncread('ADV.nc', 'dt');
tmpP = ncread('ADV.nc', 'pressure');
tmpDTsen = ncread('ADV.nc', 'dtHPR');
tmpRoll = ncread('ADV.nc', 'roll');
tmpSSC35 = ncread('ADV.nc', 'ssc35cm');
figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP, ...
             stADV.DATA.dt(iStart:iEnd), stADV.DATA.seapress(iStart:iEnd))
figure; plot(datetime(tmpDTsen,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpRoll, ...
            stVECT.dtSEN, stADV.QC.roll)
figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpSSC35, ...
            stADV.DATA.dt(iStart:iEnd), stVECT.BT2.OBS2c(iStart:iEnd))

% create global attributes (sensor description)
ncwriteatt('ADV.nc','/','description', 'Time series of (continuous) pressure, velocity, and suspended sediment concentration recorded by an acoustic doppler velocimeter co-located with two optical backscatter sensors in the surf zone at Follets Island during Hurricane Harvey.');
ncwriteatt('ADV.nc','/','comments', 'See ADV.pdf for depiction of instrument deployment setup and storm deposits.');
ncwriteatt('ADV.nc','/','creation_date', datestr(now));
ncwriteatt('ADV.nc','/','modification_date', 'none');
ncwriteatt('ADV.nc','/','institution', 'Texas A&M University, Department of Ocean Engineering');
ncwriteatt('ADV.nc','/','host', 'DesignSafe-CI');
ncwriteatt('ADV.nc','/','contact', 'kanarde@tamu.edu or figlusj@tamu.edu');
ncwriteatt('ADV.nc','/','sensor', 'Nortek Vector VEC12420 (fixed stem)');
ncwriteatt('ADV.nc','/','ADV_samp_frequency', '16 Hz');
ncwriteatt('ADV.nc','/','ADV_latitude', '29.023007801098682');
ncwriteatt('ADV.nc','/','ADV_longitude', '-95.19097459203216');
ncwriteatt('ADV.nc','/','bed_elevation_prestorm', '-0.821 m NAVD88');
ncwriteatt('ADV.nc','/','bed_elevation_poststorm', '-0.746 NAVD88');
ncwriteatt('ADV.nc','/','ADV_measurement_volume_prestorm', '-0.471 m NAVD88 (35 cm above bed)');
ncwriteatt('ADV.nc','/','ADV_PT_elevation_prestorm', '-0.046 m NAVD88 (70 cm above bed)');
ncwriteatt('ADV.nc','/','ADV_nominal_vel_range', '2.00 m/s');
ncwriteatt('ADV.nc','/','ADV_sampling_volume', '14.9 mm');
ncdisp('ADV.nc')

%% create netcdf file for PT-4

% % create variables
% iStart = 1; iEnd = 7644960; % this is the windowed data set
% nccreate('PT4.nc','pressure', 'Dimensions', {'time', length(stRRU.PT(4).seapress(iStart:iEnd))})
% nccreate('PT4.nc','dt', 'Dimensions', {'time', length(stRRU.PT(4).dt(iStart:iEnd))})
% 
% % write data to file
% ncwrite('PT4.nc', 'pressure', stRRU.PT(4).seapress(iStart:iEnd))
% ncwrite('PT4.nc', 'dt', datenum(stRRU.PT(4).dt(iStart:iEnd)+hours(5))) % convert to UTC from LDT
% 
% % create local attributes (units and description)
% ncwriteatt('PT4.nc','dt', 'units', 'datenum');
% ncwriteatt('PT4.nc','dt', 'description', 'time relative to UTC (LDT + 5 hours)');
% ncwriteatt('PT4.nc','pressure', 'units', 'dbar');
% ncwriteatt('PT4.nc','pressure', 'description', 'gauge pressure, meaning pressure corrected for atmospheric/barometric pressure; barometric pressure data was measured by a subaerial mounted PT at Follets Island (see ATM.nc).');
% 
% % testing
% tmpDT = ncread('PT4.nc', 'dt');
% tmpP = ncread('PT4.nc', 'pressure');
% %figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP)
% figure; plot(datetime(tmpDT,'ConvertFrom','datenum', 'TimeZone', 'UTC'), tmpP, ...
%              stRRU.PT(4).dt, stRRU.PT(4).seapress)
% 
% % create global attributes (sensor description)
% ncwriteatt('PT4.nc','/','description', 'Time series of gauge pressure recorded by a pressure transducer housed within a shallow PVC well located in the mid-barrier at Follets Island during Hurricane Harvey.');
% ncwriteatt('PT4.nc','/','comments', 'See PT-4.pdf for depiction of instrument deployment setup.');
% ncwriteatt('PT4.nc','/','creation_date', datestr(now));
% ncwriteatt('PT4.nc','/','modification_date', 'none');
% ncwriteatt('PT4.nc','/','institution', 'Texas A&M University, Department of Ocean Engineering');
% ncwriteatt('PT4.nc','/','host', 'DesignSafe-CI');
% ncwriteatt('PT4.nc','/','contact', 'kanarde@tamu.edu or figlusj@tamu.edu');
% ncwriteatt('PT4.nc','/','sensor', 'RBR Dwave 41454');
% ncwriteatt('PT4.nc','/','PT_samp_frequency', '16 Hz');
% %ncwriteatt('PT3.nc','/','PT_latitude', '28.616407057141966');
% %ncwriteatt('PT3.nc','/','PT_longitude', '-95.94081907817059');
% %ncwriteatt('PT3.nc','/','PT_elevation_prestorm', '0.222 m NAVD88');
% %ncwriteatt('PT3.nc','/','PT_elevation_poststorm', '0.260 m NAVD88');
% %ncwriteatt('PT3.nc','/','bed_elevation_prestorm', '0.222 m NAVD88');
% %ncwriteatt('PT3.nc','/','bed_elevation_poststorm', '0.384 m NAVD88');
% ncdisp('PT4.nc')
