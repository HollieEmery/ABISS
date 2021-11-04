%% ABISS log file reader
% Created: 2018-08-24
% Last modified: 2018-10-09
% Chris McRaven
%
% Script for loading binary log files produced by the ABISS data logger.
% See abiss_shmem.h and abiss_datalog.c for details on structure of log files.

function [timearray, dataarray, dataheader] = abiss_parse_log_file(fullpath)

%% Load file as a byte array

% use previous path if run by hand
if not(exist('dirpath', 'var'))
    dirpath = '';
end

if not(exist('fullpath','var'))
    [filename, dirpath] = uigetfile([dirpath filesep '*.dat']);
    fullpath = [dirpath filesep filename];
else
    [dirpath, filename] = fileparts(fullpath);
end
%dirpath = 'C:\Users\cmcraven\Nextcloud\WHOI\Projects\ABISS\DOCUMENTS\DATA\20181007';
%filename = 'abiss_120181007_174941.dat';
fullpath
fd = fopen(fullpath);
raw = fread(fd, Inf, '*uint8', 'n');
fclose(fd);

%% Determine file type and set up data structure sizes
% All files have the same naming; contents are written in packed structure
% format.  RGA and battery voltage are fixed length.  Optodes and "arduino"
% are variable length.  See "abiss_datalog.c" for details.
filetype = sscanf(filename, 'abiss_%c');
rawsize = size(raw);
struct = {}; % 'struct' cell array and 'bytes' array  
bytes = [];  % are used in the parse loop below

switch filetype
    case {'1', '2', 'A'} % Optodes, "arduino"; variable length
        struct = { 'double' ; 'int32' ; 'char' };
        bytes = [8, 4, 0]; 
    case 'R' % RGA, fixed length
        rgalength = 2568;
        struct = { 'double' ; 'int32' ; 'char' };
        bytes = [8, 4, rgalength];
    case 'V' % Lander battery voltage, fixed length
        struct = { 'double' ; 'single' };
        bytes = [8, 4];
    otherwise
        error([ 'Unexpected file type "' filetype '"']);
end

%% Parse the raw data
structsize = size(struct);
parsedlog = cell(1,structsize(1));
offset = 1;
strlen = 0; % used for variable-length char data
ii = 1;
while offset < rawsize(1) 
    for jj = 1:structsize(1)
        cstart = offset;
        if bytes(jj) == 0
            cend = offset + strlen - 1;
        else
            cend = offset + bytes(jj) - 1;
        end
        
        % catch when the file didn't finish writing
        try 
            sel = raw(cstart:cend);
        catch ME
            if (strcmp(ME.identifier,'MATLAB:badsubscript'))
                disp('reached end of raw data a little early, stopping raw parsing');
                offset = rawsize(1);
                break;
            end
        end
        
        % Parse data block based on structure
        switch struct{jj}
            case 'single'
                parsedlog{ii,jj} = raw2single(sel);
            case 'double'
                parsedlog{ii,jj} = raw2double(sel);
            case 'int32'
                strlen = raw2int(sel);
                parsedlog{ii,jj} = strlen;
            case 'char'
                parsedlog{ii,jj} = transpose(char(sel));
            otherwise
                error('Unexpected struct type.')
        end
        offset = cend + 1;
    end
    ii = ii + 1;
end

%% process data

switch filetype
    case 'R'
        % process the char data into an array of int32 numbers
        dataheader = {'MassSpectrum'};
        % 3rd column is the RGA data; only use when the raw length is correct 
        validentries = find(cellfun(@length,parsedlog(:,3))==rgalength);
        chararray = cell2mat(parsedlog(validentries,3));
        timearray = cell2mat(parsedlog(validentries,1));
        
        nrows = size(chararray,1);
        intarraysize = bytes(3)/4;
        dataarray = zeros(nrows, intarraysize);
        
        wbtext = sprintf('Parsing %d records of RGA data', nrows); 
        wb = waitbar(0, wbtext);
        % this is soooo slooooooowwwww
        for ii = 1:nrows
            waitbar( double(ii/nrows), wb);
            % for jj = 1:intarraysize
            %   cstart = 4*(jj-1) + 1;
            %   cend = cstart + 3;
            %   sel = transpose(chararray(ii,cstart:cend));
            %   dataarray(ii, jj) = raw2int(sel);
            % end
            charrow=chararray(ii,:);
            charrow = reshape(transpose(charrow),4,642);
            charrow = num2cell(charrow,1);
            dataarray(ii,:) = cellfun(@raw2int, charrow);
        end
        close(wb);
    case 'V'
        timearray = cell2mat(parsedlog(:,1));
        dataheader = {'BatteryVoltage[V]'};
        dataarray = cell2mat(parsedlog(:,2));
    case {'1','2'}
        timearray = cell2mat(parsedlog(:,1));
        %'MEASUREMENT %d	766	O2Concentration[uM]	16.084	AirSaturation[%]	3.927	Temperature[Deg.C]	3.988	CalPhase[Deg]	58.721	TCPhase[Deg]	59.520	C1RPh[Deg]	65.459	C2RPh[Deg]	5.939	C1Amp[mV]	1523.1	C2Amp[mV]	978.6	RawTemp[mV]	569.8??'}
        dataheader = {'OptodeModelNo','OptodeSerialNo','O2Concentration[uM]','AirSaturation[%%]','Temperature[Deg.C]','CalPhase[Deg]','TCPhase[Deg]','C1RPh[Deg]','C2RPh[Deg]','C1Amp[mV]','C2Amp[mV]','RawTemp[mV]'};
        optodestr = 'MEASUREMENT	%d	%d	O2Concentration[uM]	%f	AirSaturation[%%]	%f	Temperature[Deg.C]	%f	CalPhase[Deg]	%f	TCPhase[Deg]	%f	C1RPh[Deg]	%f	C2RPh[Deg]	%f	C1Amp[mV]	%f	C2Amp[mV]	%f	RawTemp[mV]	%f';
        %optodestr = 'MEASUREMENT %f %f O2Concentration[uM] %f AirSaturation[%] %f Temperature[Deg.C] %f CalPhase[Deg] %f TCPhase[Deg] %f C1RPh[Deg] %f C2RPh[Deg] %f C1Amp[mV] %f C2Amp[mV] %f RawTemp[mV] %f';
        nrows = size(parsedlog, 1);
        dataarray = cell(1,12);
        newtimearray = [];
        jj = 1;
        for ii = 1:nrows
            sel = parsedlog{ii,3};
            if ~isempty(sel)
                if contains(sel, 'MEASUREMENT')
                    %
                    newtimearray(jj) = timearray(ii);
                    dataarray(jj,:) = textscan(sel, optodestr);
                    jj = jj + 1;
                end
            end
        end
        timearray = transpose(newtimearray);
    otherwise
        timearray = cell2mat(parsedlog(:,1));
        dataarray = parsedlog(:,3);
        warning('Data processing not set up for this filetype. Returning raw stream.');
end

end

%% functions
function int = raw2int( raw )
hexstr = reshape(transpose(flip( dec2hex(raw, 2) )), 1, 8);
int = typecast(uint32(hex2dec(hexstr)),'int32');
end

function single = raw2single( raw )
hexstr = reshape(transpose(flip( dec2hex(raw, 2) )), 1, 8);
single = typecast(uint32(hex2dec(hexstr)),'single');
end

function double = raw2double( raw )
hexstr = reshape(transpose(flip( dec2hex(raw, 2) )), 1, 16);
double = hex2num(hexstr);
end