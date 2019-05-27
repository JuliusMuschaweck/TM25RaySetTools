function [exRot, eyRot] = RotDir( loc, radius )
    loc2 = loc(1)^2 + loc(2)^2;
    r2 = radius^2;
    if (loc2 < r2*1e-15)
        exRot = [1,0,0];
        eyRot = [0,1,0];
        return;
    end
    c = (r2-loc2)/(r2+loc2);
    mc = 1-c;
    % z = cos theta
    loc2norm = sqrt(loc2);
    fac = 1 + c;
    s = loc2norm * fac / radius;
    % s = sqrt(1-c^2);
    u = [-loc(2), loc(1)] / radius;
    u = u / norm(u);
    ux = u(1);
    uy = u(2);
    exRot = [c + ux^2*mc, ux*uy*mc, -uy*s];
    eyRot = [ux*uy*mc, c + uy^2*mc, ux*s];
end