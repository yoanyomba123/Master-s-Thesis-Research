function result = imagePreprocessing(volume)
dimensions = [5.51, 9.51, 501.98, 478.98];
    for i = 1: size(volume, 3)
        % crop image
        I = volume(:,:,i);
        I = imcrop(I, dimensions);
        % Get rid of brigh tspot
        I_t = I == 1;
        I_t = bwareaopen(I_t, 600);
        I(I_t) = 0; 
        % identify layers and get rid of B-Ground
        Im_tmp = imquantize(I, multithresh(I, 2));
        I(Im_tmp == 1) = multithresh(I, 1);
        result(:,:,i) = I;
    end
end

