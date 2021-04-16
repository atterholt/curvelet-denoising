% James Atterholt, Caltech 2020
% Curvelet median filter, originally made to remove bad traces in DAS

function MAT = MedianFilter(MAT,a,b,c)

% MAT = input matrix
% a = number of columns over which to compute the median
% b = number of rows over which to compute the median
% c = ratio over the median over which a number is considered to be an
% outlier

% Starting point: a = 50, b = 5, c = 10

absMAT = abs(MAT);

% Get a medians matrix
PartialMedians = movmedian(absMAT,a,2);
medians = movmedian(PartialMedians,b,1);

% Get a comparisons matrix
compMAT = (absMAT./medians);

% Find the bad values
[row,col] = find(compMAT > c);

% Set the bad values to their neighbors

% Note: all these if then statements exist for 2 reasons
% 1. To account for lines at the edge of the array (which only have one
% adjacent station
% 2. To account for cases where two bad traces are next to each other (a
% common occurrence)
for x = 1:length(row)
    if col(x) == 1 || col(x) == 2
        if compMAT(row(x),col(x)+1) <= c
            MAT(row(x),col(x)) = MAT(row(x),col(x)+1);
        else
            MAT(row(x),col(x)) = MAT(row(x),col(x)+2);
        end
    elseif col(x) == size(MAT,2) || col(x) == size(MAT,2)-1
        if compMAT(row(x),col(x)-1) <= c
            MAT(row(x),col(x)) = MAT(row(x),col(x)-1);
        else
            MAT(row(x),col(x)) = MAT(row(x),col(x)-2);
        end
    else
        if (compMAT(row(x),col(x)-1) <= c) && (compMAT(row(x),col(x)+1) <= 5)
            MAT(row(x),col(x)) = (MAT(row(x),col(x)-1)+MAT(row(x),col(x)+1))/2;
        elseif (compMAT(row(x),col(x)-1) <= c)
            MAT(row(x),col(x)) = (MAT(row(x),col(x)-1)+MAT(row(x),col(x)+2))/2;
        elseif (compMAT(row(x),col(x)+1) <= c)
            MAT(row(x),col(x)) = (MAT(row(x),col(x)-2)+MAT(row(x),col(x)+1))/2;
        else
            MAT(row(x),col(x)) = (MAT(row(x),col(x)-2)+MAT(row(x),col(x)+2))/2;
        end
    end
end

end