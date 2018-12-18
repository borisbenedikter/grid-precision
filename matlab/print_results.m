function [] = print_results(MAT, filename)
%print_results - Print the integration results.
%
% Syntax: [] = print_results(MAT, filename)
%

[rows, cols] = size(MAT);

fileID = fopen(filename, 'w');
fprintf(fileID, '%20s %20s %20s %20s %20s %20s %20s \n', ...
    'Time [s]', 'X [km]', 'Y [km]', 'Z [km]', 'VX [km]', 'VY [km]', 'VZ [km]');
for j = 1:rows
    fprintf(fileID, '%20.9e %20.9e %20.9e %20.9e %20.9e %20.9e %20.9e \n', MAT(j, :));
end

end
    