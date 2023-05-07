classdef LuminanceInterpolator < handle

    properties
        fac = 1
        boundingbox = struct()
        gridsize = struct()
        gridx = []
        gridy = []
        gridkx = []
        gridky = []
        data = []
    end

    methods

        function Read4DTable(obj, filename)
            fh = fopen(filename,'r');
            
            obj.boundingbox.lo = fread(fh,4,'float32');
            obj.boundingbox.hi = fread(fh,4,'float32');
            
            obj.gridsize.nx = fread(fh,1,'uint64');
            obj.gridsize.ny = fread(fh,1,'uint64');
            obj.gridsize.nkx = fread(fh,1,'uint64');
            obj.gridsize.nky = fread(fh,1,'uint64');
            obj.gridx = linspace( obj.boundingbox.lo(1), obj.boundingbox.hi(1), obj.gridsize.nx);
            obj.gridy = linspace( obj.boundingbox.lo(2), obj.boundingbox.hi(2), obj.gridsize.ny);
            obj.gridkx = linspace( obj.boundingbox.lo(3), obj.boundingbox.hi(3), obj.gridsize.nkx);
            obj.gridky = linspace( obj.boundingbox.lo(4), obj.boundingbox.hi(4), obj.gridsize.nky);

            np = obj.gridsize.nx * obj.gridsize.ny * obj.gridsize.nkx * obj.gridsize.nky;
            raw = fread(fh,np,'float32');
            obj.data = nan(obj.gridsize.nx, obj.gridsize.ny, obj.gridsize.nkx, obj.gridsize.nky);
            idx = 0;
            for ix = 1:obj.gridsize.nx
                for iy = 1:obj.gridsize.ny
                    for ikx = 1:obj.gridsize.nkx
                        for iky = 1:obj.gridsize.nky
                            idx = idx + 1;
                            obj.data(ix,iy,ikx,iky) = raw(idx);
                        end
                    end
                end
            end
            fclose(fh);
        end

        function rv = BoundingBox(obj)
            rv = obj.boundingbox;% return struct with xmin, xmax, ymin, ymax
        end

        function rv = Scale(obj, fac)
            % return previous value
            rv = obj.fac;
            obj.fac = fac;
        end

        function rv = IsInBoundingBox(obj, x,y,kx,ky)
            if x < obj.boundingbox.lo(1) || x > obj.boundingbox.hi(1) ...
                 || y  < obj.boundingbox.lo(2) || y  > obj.boundingbox.hi(2) ...
                 || kx < obj.boundingbox.lo(3) || kx > obj.boundingbox.hi(3) ...
                 || ky < obj.boundingbox.lo(4) || ky > obj.boundingbox.hi(4)
                rv = false;
            else
                rv = true;
            end
        end

        function rv = Luminance(obj,  x, y, kx, ky) % x, y, kx, ky must be vectors of same orientation and length
            tmp = interpn(obj.gridx, obj.gridy, obj.gridkx, obj.gridky, obj.data, x, y, kx, ky,'linear',0);
            rv = tmp * obj.fac;
        end

    end

end