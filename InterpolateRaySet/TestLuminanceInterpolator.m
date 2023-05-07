clear;
lInt = LuminanceInterpolator;
lInt.Read4DTable("LumDiskTest_luminanceTable.dat");

test1 = lInt.Luminance(0,0,0,0);

ix = -5:0.025:5;
iy = -5:0.025:5;

ikx = 0;
iky = 0;

[xx,yy] = meshgrid(ix, iy);
kkx = ikx * ones(size(xx));
kky = iky * ones(size(xx));

val = lInt.Luminance(xx(:), yy(:), kkx(:), kky(:));
val = reshape(val, size(xx));

figure(1);
imagesc(ix,iy,val);
clim([0,max(val(:))+1e-7]);
colorbar;

figure(2);
val2 = lInt.data(1:lInt.gridsize.nx,1:lInt.gridsize.ny,round((lInt.gridsize.nkx)/2),round(lInt.gridsize.nky/2));
imagesc(val2);
colorbar;

tmp = lInt.data(:);
%%
min(tmp)
max(tmp)
mean(tmp)
figure(3);
histogram(tmp,2000);

sum(tmp ~= 0) / length(tmp)

stmp = sort(tmp);
figure(4);
plot(log10(stmp+1e-3));