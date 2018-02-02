function [ spectrum wavelength ] = getOsaTrace()
%GETOSATRACE gets wavelength and spectrum data from the ando aq6317 OSA
%trace B

% Find a GPIB object.
OSA = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 5, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(OSA)
    OSA = gpib('NI', 0, 5);
else
    fclose(OSA);
    OSA = OSA(1);
end

% set input buffer size so MATLAB can read the (long) char arrays output
% from OSA
set(OSA, 'InputBufferSize', 80192);

% Connect to instrument object, OSA.
fopen(OSA);

% Get spectrum data (change B to A or C for other traces)
spectrum = query(OSA, 'LDATB');

% flushinput(OSA);

% Get wavelength data
wavelength = query(OSA, 'WDATB');

% output data is char arrays; clean it up and convert to doubles

spectrum = str2num(spectrum(6:end));
wavelength = str2num(wavelength);

spectrum = spectrum(2:end-1);
wavelength = wavelength(3:end-1);

% Disconnect from instrument object, obj1.
fclose(OSA);


end

