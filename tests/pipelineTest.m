clearvars;

javaaddpath('../src/vision7_symphony/Vision.jar');

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import edu.ucsc.neurobiology.vision.util.*;

% Try to parse the file
% fileName = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220329C_mouse/data052/data052000.bin'; %'data000000.bin';

if ispc
    fileName = 'C:\Users\Public\Documents\Data\data000000.bin';
else
    fileName = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220406C/data011/data011000.bin';
end

streamPort = 9000; % 9876; 9000

fileID = fopen(fileName);
% Read the whole file.
A = fread(fileID);
fclose(fileID);

% f = java.io.FileInputStream('data000000.bin');
%%

host = InetAddress.getLocalHost();
socket = Socket(host, streamPort);
output = socket.getOutputStream();
% Send the spike finding command.
% output.writeInt(int32(34));
output.write( java.nio.ByteBuffer.allocate(4).putInt(34).array() );
for jj = 1 : length(A)
    output.write(A(jj));
end

output.close();
socket.close();

