function DrawBounds(I2,I3)
[B,L] = bwboundaries(I2,'noholes');
figure;
imshow(I3);
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end

end

