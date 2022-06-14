clear;
fileID = fopen('characteristicCurve.dat');
%format: 8 byte int for nCells, 2 by n matrix of 4 byte floats, row major, first row is etendue, second row is luminance.
nCells = fread(fileID,1,'uint64');
etendue = fread(fileID,nCells,'float');
luminance = fread(fileID,nCells,'float');
fclose(fileID);

%%
figure(1);
clf;
plot(etendue,luminance);