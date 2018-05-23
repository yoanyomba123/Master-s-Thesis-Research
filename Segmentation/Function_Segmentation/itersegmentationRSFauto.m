function [energyField, stopNormSeg] = itersegmentationRSFauto(imgInput, Parameter)

Img = imgInput;
c0 = 2;
iterNum = Parameter.iterNum;
lambda1 = Parameter.lambda1;
lambda2 = Parameter.lambda2;
nu = Parameter.nu*255*255; % coefficient of the length term
initialLSF = Parameter.initialLSF;
% logicLSF = roipoly(imgInput,initialLSFCoor(1,:),initialLSFCoor(2,:));
initialLSF = initialLSF .* c0;
% initialLSF = initialLSF .* (~logicLSF) + ones([512 512]).* -c0 .* logicLSF;
figureBose = Parameter.figureBose;

u = initialLSF;
stopNormSeg = ones(1, iterNum);
if (figureBose == 1)
	figure;imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
	title('Initial contour');
	[~,~] = contour(u,[0 0],'r');
	pause(0.1);
end

timestep = 0.1;% time step
mu = 1;% coefficient of the level set (distance) regularization term P(\phi)

epsilon = 1.0;  % the papramater in the definition of smoothed Dirac function
sigma = Parameter.sigma;    % scale parameter in Gaussian kernel
% Larger, more robust
% Smaller, more accurate
K=fspecial('gaussian',round(2*sigma)*2+1,sigma);     % the Gaussian kernel
I = Img;
KI=conv2(Img,K,'same');     % compute the convolution of the image with the Gaussian kernel outside the iteration
                            % See Section IV-A in the above IEEE TIP paper for implementation.
                                                 
KONE=conv2(ones(size(Img)),K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
                                       % See Section IV-A in the above IEEE TIP paper for implementation.

% start level set evolution
for n=1:iterNum
    uPrevious = u;
    u=RSF(u,I,K,KI,KONE, nu,timestep,mu,lambda1,lambda2,epsilon,1);
    stopNormSeg(1, n) = norm(double(uPrevious == u), 2);
    if(stopNormSeg(1, n) <= 2); break; end
    if mod(n,20)==0
    	if (figureBose == 1)
	        pause(0.01);
	        imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
	        [~,~] = contour(u,[0 0],'r');
	        iterNum=[num2str(n), ' iterations'];
	        title(iterNum);
	        hold off;
	    end
    end
end
if (figureBose == 1)
	imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
	[~,~] = contour(u,[0 0],'r');
	totalIterNum=[num2str(n), ' iterations'];
	title(['Final contour, ', totalIterNum]);
	% figure;
	% mesh(u);
	% title('Final level set function');
end
energyField = u;
end



