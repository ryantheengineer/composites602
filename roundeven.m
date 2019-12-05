function [outval] = roundeven(inval)
% Round inval to the nearest even integer
    if inval < 2
        error('roundeven function only takes values above 2');
    end
    
    extra = mod(inval,2);
    
    if extra > 1
        outval = inval + (1-(extra-1));
    else
        outval = inval - extra;
    end

end