function [ var ] = LoadParameter( file_in )
% Load mesh generation parameter
    [tmps,var]=textread(file_in,'%s %f',5);
end

