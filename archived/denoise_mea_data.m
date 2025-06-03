
Fs = 2e4;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);

60*[3,5,7,9,11,13,15]
           
data_path = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220607C/data000';

% Open the Litke/Symphony data folder for reading.
h = edu.ucsc.neurobiology.vision.io.RawDataFile(data_path);

% Get the header.
header = h.getHeader();

channel = 1;
data = h.getData(channel, 0, header.getNumberOfSamples());

% Pass through the butterworth filter to remove 60 cycle noise.
dfilt = filtfilt(d, double(data));

[popen,fopen] = periodogram(double(data),[],[],Fs);
[pbutt,fbutt] = periodogram(dfilt,[],[],Fs);

plot(fopen,20*log10(abs(popen)),fbutt,20*log10(abs(pbutt)),'--')