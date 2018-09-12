function [] = setWriteToFile(arg)
%This specifies whether you want to write your window outputs to a text
%file 
% arg = 0, means no
% arg = 1, means yes
writeToFile = arg; 
save writeToFile
end

