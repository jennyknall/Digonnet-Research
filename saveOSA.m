%% Instrument Connection

% Find a GPIB object.
osa = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', 5, 'Tag', '');

% Create the GPIB object if it does not exist
% otherwise use the object that was found.
if isempty(osa)
    osa = gpib('NI', 0, 5);
else
    fclose(osa);
    osa = osa(1);
end

% set input buffer size and timeout 
set(osa, 'InputBufferSize', 100000);
set(osa, 'Timeout', 1000.0);

% Connect to instrument object, osa
fopen(osa);

fprintf(osa, 'RESln 0.01') 
resoll=eval(query(osa, 'RESln?'));

fprintf(osa, 'SRQ1');
[a,statusByte]=spoll(osa); % change if your trace is not trace A - to keep it simple, set your active trace from the front panel of the OSA to Trace A. Block the other two traces B and C. 
fprintf(osa,'SGL');

statusByte=0;
disp('waiting')

while(mod(statusByte,2)==0)
    [a,statusByte]=spoll(osa);
end

disp('sweep done')
% getting the power data
power1=query(osa,'LDATA');
power1=str2num(power1);
power1=power1(2:end);
% getting the wavelength data
wave1=query(osa,'WDATA');
wave1=str2num(wave1);
wave1=wave1(2:end).*1e-9;
disp('done reading')

figure(1)
plot(wave1,power1)
Hold on


pathroot='F:\My_Files\Mina_data\data_10_26_2017'; %customize it for your path 
FName=['1020nm_fiber_laser_spec_1mW_pout']; % customize it for your own labeling 
save([ pathroot '\SPECTRA_' FName '_' int2str(fix(clock))], 'wave1', 'power1', 'resoll')

fcolse(osa);
