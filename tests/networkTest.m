clearvars;

addpath('../src/simulation/');
javaaddpath('../src/vision7_symphony/Vision.jar');

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import edu.ucsc.neurobiology.vision.util.*;

fileID = fopen('data000000.bin');
% Read the whole file.
A = fread(fileID);
fclose(fileID);

% f = java.io.FileInputStream('data000000.bin');

host = InetAddress.getLocalHost();
socket = Socket(host, 9001);
output = socket.getOutputStream();
for jj = 1 : length(A)
    output.write(A(jj));
end

output.close();
socket.close();


return

% Create a SymphonyStream
% n = edu.ucsc.neurobiology.vision.io.SymphonyStream();

fileID = fopen('data000000.bin');

% Read the whole file.
A = fread(fileID);

% onebyte = fread(fileID,4,'*ubit8');

fclose(fileID);

f = java.io.FileInputStream('data000000.bin');

h = edu.ucsc.neurobiology.vision.io.RawDataHeader512(f);

fileName = [char(h.getExperimentIdentifier()), '\',char(h.getDatasetName()), '.bin'];


%%
HEADER_LENGTH_TAG = 0;
TIME_TAG = 1;
COMMENT_TAG = 2;
FORMAT_TAG = 3;
ARRAY_ID_TAG = 4;
FREQUENCY_TAG = 5;
TRIGGER_TAG = 6; %//deprecated on 7/3/2015, mgrivich
DATASET_IDENTIFIER_TAG = 7;
TRIGGER_TAG_V2 = 8; %//new on 7/3/2015
DATA_TAG = 499;
FILE_TYPE = 0x512;

fileID = fopen('data000000.bin');

tagsLeft = true;
while (tagsLeft)
    tag = fread(fileID,1,'int');
    len = fread(fileID,1,'int');
    switch tag
        case HEADER_LENGTH_TAG
            headerLength = fread(fileID,1,'int');
        case TIME_TAG
            timeBase = fread(fileID,1,'int');
            secondsTime = fread(fileID,1,'long');
        case COMMENT_TAG
%             b = new byte[len];
%             for (int i = 0; i < length; i++) b[i] = stream.readByte();
%             this.comment = new String(b);
%             if (DEBUG) System.out.println("comment: " + comment);
%             break;
        case FORMAT_TAG
            this.format = fread(fileID,1,'int');

        case ARRAY_ID_TAG
            this.nElectrodes = fread(fileID,1,'int');
            this.arrayID = fread(fileID,1,'int');

        case FREQUENCY_TAG
            this.frequency = fread(fileID,1,'int');

        case TRIGGER_TAG
            fread(fileID,1,'int');
            fread(fileID,1,'int');

        case TRIGGER_TAG_V2
            fread(fileID,1,'int');
            fread(fileID,1,'int');
            fread(fileID,1,'int');
            fread(fileID,1,'int');

        case DATASET_IDENTIFIER_TAG
%             b = new byte[length];
%             for (int i = 0; i < length; i++) b[i] = stream.readByte();
%             this.datasetIdentifier = new String(b);
%             if (DEBUG) System.out.println("datasetIdentifier: " + datasetIdentifier);
%             break;

        case DATA_TAG
            this.nSamples = fread(fileID,1,'int');
%             if (DEBUG) System.out.println("nSamples: " + nSamples);
%             tagsLeft = false; // DATA_TAG signals the end of the header
%             break;
        otherwise
%             b = new byte[length];
%             for (int i = 0; i < length; i++) b[i] = stream.readByte();
%             this.comment = new String(b);
%             System.out.println("Warning: Unknown Tag:" + tag + " Length:" + length);
%             break;
    end
end

% onebyte = fread(fileID,4,'*ubit8');

fclose(fileID);


return
javaaddpath('./Vision.jar')


n = edu.ucsc.neurobiology.vision.io.NetworkDaemon();


