function Volume = Convert_to3d(cell_array)

for i = 1: length(cell_array)
   Volume(:,:,i) = gather(cell_array{i}); 
end
end

