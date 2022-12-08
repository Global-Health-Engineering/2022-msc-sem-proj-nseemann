tic;
disp('timer started');

%% Distance matrix
fileID = fopen('../../data/interm_data/dist_matrix_ID_filtered.csv');
textscan(fileID,'%s',3,'Delimiter',',');
C_dist = textscan(fileID,'%d %d %f','Delimiter',',');
fclose(fileID);


%Convert into a matrix
dist_mat = zeros(length(max(max(C_dist{1},max(C_dist{2})))));
for i = 1:length(C_dist{1})
    dist_mat(C_dist{1}(i),C_dist{2}(i)) = C_dist{3}(i);
end

%% 