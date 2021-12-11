function [modified_array] = reshapearray(cruncharray)

if length(cruncharray(1,:)) > 1
    modified_array = cruncharray';
else
    modified_array = cruncharray;
end