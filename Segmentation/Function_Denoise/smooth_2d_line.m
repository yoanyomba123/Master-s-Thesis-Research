function img = smooth_2d_line(img, power)
    % SG function implemented line by line by the system
    for i = 1:1:size(img, 1)
        img(i, :) = smooth(img(i, :), power, 'sgolay');
    end

    for i = 1:1:size(img, 1)
        img(:,i) = smooth(img(:,i), power, 'sgolay');
    end
end