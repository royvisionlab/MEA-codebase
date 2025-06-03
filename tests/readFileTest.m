clearvars;

javaaddpath('../src/vision7_symphony/Vision.jar');

% import java.io.*;
% import java.net.*;
% import java.util.ArrayList;
% import java.util.Collections;
% import java.util.List;
import edu.ucsc.neurobiology.vision.util.*;

% Try to parse the file
% fileName = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220329C_mouse/data052/data052000.bin'; %'data000000.bin';
% fileName = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220406C/data011/data011000.bin';
fileName = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220406C/data033/data033000.bin';

h = edu.ucsc.neurobiology.vision.io.RawDataFile(fileName);

header = h.getHeader();

startSample = 0;
% foo = h.getData(startSample, header.getNumberOfSamples());
foo = h.getData(startSample, min(400000,header.getNumberOfSamples()));

plot(foo(:,162))





