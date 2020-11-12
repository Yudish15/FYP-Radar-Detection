function data = read_BRDF_binary(binary_file,num_r_bins,pulses_per_file)

% Function to read .brd files

fileID = fopen(binary_file);
data = single(fread(fileID,Inf,'single'));
fclose(fileID);

data = reshape(data(1:2:end) + 1i.*data(2:2:end),[num_r_bins,pulses_per_file]);

