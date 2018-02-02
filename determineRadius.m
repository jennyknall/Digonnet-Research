function [radius] = determineRadius(loops,r_clad)

switch loops
    case 1
        radius = r_clad;
    case 2
        radius = 2*r_clad;
    case 3
        radius = (1+2/sqrt(3))*r_clad;
    case 4
        radius = (1+sqrt(2))*r_clad;
    case 6
        radius = (1+sqrt(7))*r_clad;
    case 7
        radius = 3*r_clad;
    case 9
        radius = (1+sqrt(12))*r_clad;
    case 13
        radius = (1+sqrt(12))*r_clad;
    otherwise
        radius = nan;
end

end