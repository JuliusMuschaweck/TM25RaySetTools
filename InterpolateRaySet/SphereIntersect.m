function [p, lambda] = SphereIntersect(center, radius, r0, k)
    r0c = r0 -center;
    r0c = r0c(:);
    k = k(:);
    a = k' * k;
    b = k' * r0c; 
    c = r0c'*r0c - radius^2;
    D = b^2 - a * c;
    if (D>=0)
        lambda = 1/(a) * (-b + sqrt(D));
        p = r0(:) + lambda * k;
    else
        lambda = NaN;
        p = [NaN, NaN, NaN];
    end
end