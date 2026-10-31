% According to Holger Trahe, Lumileds uses the wrong NaN in at least some
% of their TM25 ray files, in the flux field.
% He said that instead of the correct 7F 80 00 01 hex for NaN, they reverse
% the byte order, as if it were big endian. That results in 
% 01 00 80 7F, which is denormalized 2.3602e-38 in little endian.
% This code reads the file which I created with Visual Studio binary editor
% to test it. Holger is right about what the byte sequences mean.

fn = "01234567_NaN_denorm.bin";
fh = fopen(fn,"r");
head = fread(fh,8,"char*1");
mydenorm = fread(fh,1,"single")
myNaN = fread(fh,1,"single")
fclose(fh);
