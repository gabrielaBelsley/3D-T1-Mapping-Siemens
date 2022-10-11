function retval = extractFromTok(x)

    if isempty(x)
        retval = [];
    else
        for i=1:length(x)
            retval(i) = str2num(x{i}{1});
        end
    end
