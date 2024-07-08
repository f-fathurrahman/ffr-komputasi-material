%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mem_str, nXB, XB] = my_print_mem(nBytes)
% given memory size in Bytes, print in appropriate units

if nBytes < 0
    error('Memory size must be non-negative!');
end
scale = floor(log(nBytes)/log(1024));
if nBytes == 0, scale = 0; end
switch scale
    case 0
        nXB = nBytes; XB = 'B';
    case 1
        nXB = nBytes/(1024); XB = 'kB';
    case 2
        nXB = nBytes/(1024^2); XB = 'MB';
    case 3
        nXB = nBytes/(1024^3); XB = 'GB';
    case 4
        nXB = nBytes/(1024^4); XB = 'TB';
    case 5
        nXB = nBytes/(1024^5); XB = 'PB';
    otherwise
        % use PB for all larger mem size
        nXB = nBytes/(1024^5); XB = 'PB';
end

if scale == 0
    mem_str = sprintf('%7.0f %s',nXB,XB);
else
    mem_str = sprintf('%7.2f %s',nXB,XB);
end

end