clear;
fileID = fopen('skewness_z_Curve.dat');
%format: 4 byte floats, 3 floats for point, 3 floats for direction,  4 byte int for nBins
% nBins+1 (!!) floats for skewness, nBins floats for dU_ds, nBins floats for dPhi_ds
axis_p = fread(fileID,3,'float');
axis_dir = fread(fileID,3,'float');
nBins = fread(fileID,1,'uint32');
skewness = fread(fileID,nBins+1,'float');
dU_ds = fread(fileID,nBins,'float');
dPhi_ds = fread(fileID,nBins,'float');
fclose(fileID);
%%
[ss,uu] = barData(skewness, dU_ds);
figure(1);
clf;
plot(ss,uu);
title('Skewness distribution of etendue');
xlabel('skewness');
ylabel('dU/ds');

%%
[ss,phiphi] = barData(skewness, dPhi_ds);
figure(2);
clf;
plot(ss,phiphi);
title('Skewness distribution of flux');
xlabel('skewness');
ylabel('d\Phi/ds');


%%
function [xx,yy] = barData(x,y)
    nx = length(x);
    ny = length(y);
    if nx ~= (ny+1)
        error('barData: wrong size');
    end
    xx = zeros(1,2*nx);
    yy = xx;
    for i = 1:nx
        xx(2*i-1) = x(i);
        xx(2*i) = x(i);
    end
    yy(1) = 0;
    for i = 1:ny
        yy(2*i) = y(i);
        yy(2*i+1) = y(i);
    end
    yy(2*nx) = 0;
end