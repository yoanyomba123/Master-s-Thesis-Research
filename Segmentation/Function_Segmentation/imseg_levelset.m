function [M1, M2, M3] = imseg_levelset(img, ParaSeg, booDebug)
    % Segmentation method based on level set.
    % Originally brought up by Chuming Li.
    % Modified by Pengwei Wu, under the supervision of Tianye Niu.

    % nu is used to  smooth the curve.
    % mu has a similar effect.
    % sigma is like the radius of the neighborhood.
    % larger sigma means that the function becomes closer to the original level set function.




    booBiasCorrection = 0;
    % Only open it when you actually need bias correction.

    A = 255; % (well, for double image, 255 pixel intervals).
    % Note that this parameter can not be changed.

    iter_outer = ParaSeg.iterOuter; % times of iteration for outer loop
    sigma = ParaSeg.sigma; % scale parameter
    timeStep = ParaSeg.timeStep;
    mu = ParaSeg.muBase / timeStep; % (u in the passage)

    nu = ParaSeg.nuBase * A ^ 2; % weight of length term (v) % 0.001 at first
    % c0 = 1;
    epsilon = ParaSeg.epsilon;
    scale = ParaSeg.scale;

    img = imresize(img, scale);

    if(size(img, 1) > 256 || size(img, 2) > 256)
        warning(['You should better down sample the image']);
    end

    % Now we began our actual segmentation
    Mask = (img > 5);





    % [~, ~] = size(img);
    % numframe = 0;

    figure;
    % the window will always be [0 255]
    imagesc(img, [0 255]); colormap(gray); hold on; axis off; axis equal;
    title('Segmentation begins'); maxscreen();

    % initialization of bias field and level set function
    b = ones(size(img));
    initialLSF(:,:,1) = randn(size(img)); % randomly initialize the level set functions
    initialLSF(:,:,2) = randn(size(img)); % randomly initialize the level set functions
    initialLSF(:,:,1) = Mask; % remove the background outside the mask
    u = sign(initialLSF);

    [~,~] = contour(u(:,:,1),[0 0],'r');
    [~,~] = contour(u(:,:,2),[0 0],'b');

    hold off
    % Kernel for clustering
    % we have implemented two different kernels
    % Will, you can turn off the kernel by the boolean followed.

    booGaussianKernel = 1;
    if(booGaussianKernel)
        Ksigma = fspecial('gaussian', round(2 * sigma) * 2 + 1, sigma); % Gaussian kernel
    else
        disk_radius = 7;
        Ksigma=fspecial('disk',disk_radius); % an alternative kernel as a truncated uniform function
    end
    KONE = conv2(ones(size(img)), Ksigma, 'same');

    totaltime = 0;
    disp(['Total time is now : ', num2str(totaltime)]);
    for n = 1:iter_outer

        t0 = cputime;
        [u, b, C] =  lse_bfe_3Phase(u, img, b, Ksigma, KONE, nu, timeStep, mu, epsilon);
        t1 = cputime;
        totaltime = totaltime + t1 - t0;

        if(mod(n,3) == 0)
            pause(0.01);
            imagesc(img,[0 255]); colormap(gray); hold on; axis off; axis equal;
            [~,~] = contour(u(:,:,1), [0 0], 'r');
            [~,~] = contour(u(:,:,2), [0 0], 'b');
            iterNum=[num2str(n), ' iterations'];
            title(iterNum);
            hold off;
        end

    end
    disp(['Total time is now : ', num2str(totaltime)]);

    if(~isempty(ParaSeg.save))
        saveas(gcf,ParaSeg.save)
    end

    H1 =  Heaviside(u(:,:,1), epsilon);
    H2 =  Heaviside(u(:,:,2), epsilon);
    M1 = H1 .* H2;
    M2 = H1 .* (1 - H2);
    M3 = (1 - H1);

    M1 = imresize(M1, 1 / scale);
    M2 = imresize(M2, 1 / scale);
    M3 = imresize(M3, 1 / scale);

    imgSeg = C(1) * M1 + C(2) * M2 + C(3) * M3;  % three regions are labeled with C1, C2, C3
    if(booDebug); figure; imagesc(imgSeg); axis off; axis equal; title('Segmented regions');
    colormap(gray); maxscreen(); end

    % if(booDebug);
    % figure;
    % imagesc(img, [0 255]); colormap(gray);hold on; axis off; axis equal;
    % [~, ~] = contour(u(:,:,1),[0 0],'r','LineWidth',1);
    % [c, h] = contour(u(:,:,2),[0 0],'b','LineWidth',1);
    % maxscreen();
    % end

    if(booBiasCorrection)
        if(booDebug); figure;
        imshow(uint8(Mask .* normalize01(b) * 200), [0 255]); colormap(gray); ...
            hold on; axis off; axis equal
        title('Estimated bias field (on the mask)');
        end

        if(booDebug); figure; title('Histogram');
        subplot(1,2,1)
        [N, X] = hist(img(:), 30); plot(X, N, 'b', 'linewidth', 2); title('Histogram of original image');
        subplot(1,2,2)
        [N, X] = hist(img_corrected(:), 30); plot(X, N, 'r', 'linewidth', 2); title('Histogram of bias corrected image');
        set(gcf, 'color' ,'w');
        end
    end

end