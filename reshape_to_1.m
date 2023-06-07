function [reshaped_matrix, total_sizes]=reshape_to_1(input_matrix)
%get the nx1 dimension for a matrix of any size and put it into that
%dimensionality.
%Compiled by NR for general use on 9/16/2013

total_sizes=size(input_matrix);

full_size=1;

for i =1:length(total_sizes)
    full_size=full_size*total_sizes(i);
end

reshaped_matrix=reshape(input_matrix,full_size, 1);

end