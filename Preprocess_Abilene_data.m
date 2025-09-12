% Read Abilene traffic matrices (5-minute intervals).
% Copy one week of data into abilene_dat for 2004/04/02–2004/04/08.

clear variables

% Constants
Ar = 12;          % Number of Abilene sites (matrix is 12x12)
At = 12*24*7;     % Number of time points for 7 days
%At = 12*24*10;    % Number of time points for 10 days 
%At = 12*24*11;   % Number of time points for 11 days
skip = 13;       % Number of header lines inside each Abilene CSV file

img_folder = './Abilene/2004/Measured';
addpath(img_folder);

% === Minimal change: list CSV files only ===
file_inf = dir(fullfile(img_folder, '*.dat'));

X = zeros(Ar, Ar, At);

i = 1;
for fid = 1:At
    % Print progress (file index and name)
    sprintf('%04d, %s', i, file_inf(fid).name)

    % Get file name (folder already on path)
    openfile = [file_inf(fid).name];

    % Read CSV skipping header lines
    M_temp = csvread(openfile, skip);

    % Store into the 3D traffic tensor
    X(:, :, i) = M_temp;
    i = i + 1;
end

save abilene_tm_2016.mat X;
%save abilene_tm_2880.mat X;
%save abilene_tm_3168.mat X;