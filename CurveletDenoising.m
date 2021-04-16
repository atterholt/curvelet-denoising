% Written by James Atterholt, 2021, using the Curvelab package

% This is, as its name suggests, the filtering script for Curvelet related
% Problems. We write it as a function so that we may parallelize this in
% separate MATLAB tasks.
% Here, we'll:
% 1. Enter the Curvelet domain.
% 2. Perform Gaussian denoising using soft thresholding
% 3. Perform velocity filtering using the standard FK methodology
% 4. Return from whence we came

% Choice = 0 (just stochastic) or 1 (just coherent) or 2 (both)

function FilteredMAT = CurveletDenoising(MAT,E,v,choice,nbscales,nbangles)

    % 2nd arg: 1 = real valued curvelets for speed
    % 3rd arg: 2 = wavelets at finest scale for speed
    % We use the wrapping routine instead of USFFT for speed.
    % Everything is for speed!
    
    C = fdct_wrapping(MAT,1,2,nbscales,nbangles);
    
    if (choice == 0) || (choice == 2)

        for s = 2:(length(C))
        
            for w = 1:length(C{s})
            
                % First do a hard threshold
                C{s}{w} = C{s}{w}.* (abs(C{s}{w}) > (E{s}{w}));
            
                % Then soften the existing coefficients
                C{s}{w} = sign(C{s}{w}).*abs(abs(C{s}{w})-abs(E{s}{w}));
            
            end
        
        end
        
    end
    
    % Now that the thresholding is over, we apply velocity filtering for
    % the velocity range that you've chosen. We do a search for each wedge
    % at each seismological scale. You will need to edit this part too if
    % you are changing the polar tiling structure.
    
    if (choice == 1) || (choice == 2)
    
        % Start your loop my friend
        for s=1:length(C)
            for w = 1:length(C{s})
            
                % For the west traveling component 
                if (s == 3 && ismember(w,v))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 4 && ismember(w,v))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 5 && ismember(w,v*2-1))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 5 && ismember(w,v*2))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 6 && ismember(w,v*2-1))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 6 && ismember(w,v*2))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 7 && ismember(w,v*4-3))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 7 && ismember(w,v*4-2))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 7 && ismember(w,v*4-1))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));
                
                elseif (s == 7 && ismember(w,v*4))
                    C{s}{w}(:) = zeros(size(C{s}{w}(:),1),size(C{s}{w}(:),2));

                % For all else
                else
                    C{s}{w}(:) = C{s}{w}(:);
    
                end
            
            end
        end
        
    end
    
    % Perform the inverse curvelet transform 
    FilteredMAT = real(ifdct_wrapping(C,1));
    
    end