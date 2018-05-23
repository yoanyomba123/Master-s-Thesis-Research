function inQueue = grow_queue(img, queue, simiTh, forbid, pixValTh)
    % img : inputted image
    % queue : initial points in queue
    %       eg : [1 3; 4 5] (x(1) y(1); x(2) y(2))
    %       row first, then column 
    % simiTh : similarity threshold
    % forbid : forbidden region
    % pixValTh : pixel value threshold
    
    % inQueue : the outputted image ('1' for the points in queue)
    % Pengwei Wu
    
    if(nargin < 5)        
        pixValTh = min(img(:));    
    end
    if(nargin < 4)
        forbid = zeros(size(img));      
    end
    if(nargin < 3)
        error('Not enough arguments in growing');
    end

    % forbid is the forbidden region of growing method
    [m, n] = size(img);

    inQueue = zeros(size(img));
    inQueue(sub2ind(size(img), queue(:,1), queue(:,2))) = 1; % all the boundary is now in queue for further growing

    while (~isempty(queue))
        t = queue(1,:);
        x = t(1); y = t(2);
        for i = max(1, x-1):min(m, x+1) % in general -- x-1:x+1
            for j = max(1, y-1):min(n,y+1)
                % cross growing
%                 if( (i==x)&&(j==y+1) || (i==x)&&(j==y-1) || (i==x-1)&&(j==y) || (i==x+1)&&(j==y) && strcmp(type,'cross') )
                    if ((abs(img(i, j) - img(x, y)) < simiTh)  && forbid(i, j) == 0 && img(i, j) >= pixValTh)
                        new = [i j];
                        if inQueue(i, j) == 0 % not a member for now
                            queue(size(queue, 1)+1, :) = new;
                            inQueue(i, j) = 1;
                        end
                    end
%                 end
            end
        end
        queue(1,:) = [];
    end  
        
end