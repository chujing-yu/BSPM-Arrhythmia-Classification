clear 
clc

% Labels=[1 2 2 3 3 3 3 ...
%         0 0 1 2 3 2 2 ...
%         0 2 3 3 2 3 2 ...
%         1 2 2 2 2 3 3 ...
%         0 3 1 3 0 3 3 ...
%         3 0 3 3 0 0 0 ...
%         3 0 0 3 0 0 3 ...
%         3 3 3 2 0 3 0 ...
%         0 2 3 2 2 1 3 ...
%         3 0 0 2 2 0 0 ...
%         3 2 0 2 2 2 2 ...
%         3 0 3 3 3 3 ...
%         2 0 2 3 0 2 2 ...
%         2 0 3 2 0 0 0 ...
%         3 0 0 3 0]

Labels=[0 1 1 2 2 2 2 ...            //0:else 1:pAF 2:cAF
        0 0 0 1 2 1 1 ...
        0 1 2 2 1 2 1 ...
        0 1 1 1 1 2 2 ...
        0 2 0 2 0 2 2 ...
        2 0 2 2 0 0 0 ...
        2 0 0 2 0 0 2 ...
        2 2 2 1 0 2 0 ...
        0 1 2 1 1 0 2 ...
        2 0 0 1 1 0 0 ...
        2 1 0 1 1 1 1 ...
        2 0 2 2 2 2 ...
        1 0 1 2 0 1 1 ...
        1 0 2 1 0 0 0 ...
        2 0 0 2 0]
% Labels=[0 1 1 1 1 1 1 ...
%         0 0 0 1 1 1 1 ...
%         0 1 1 1 1 1 1 ...
%         0 1 1 1 1 1 1 ...
%         0 1 0 1 0 1 1 ...
%         1 0 1 1 0 0 0 ...
%         1 0 0 1 0 0 1 ...
%         1 1 1 1 0 1 0 ...
%         0 1 1 1 1 0 1 ...
%         1 0 0 1 1 0 0 ...
%         1 1 0 1 1 1 1 ...
%         1 0 1 1 1 1 ...
%         1 0 1 1 0 1 1 ...
%         1 0 1 1 0 0 0 ...
%         1 0 0 1 0]
save('Labels.mat','Labels','-mat');