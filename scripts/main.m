if exist('main.m','file')
    basedir='..';
    datadir='..';
elseif exist('/code/scripts/main.m','file') %code ocean
    basedir='/code';
    datadir='/';
elseif exist('scripts','dir')
    basedir='.';
    datadir='.';
else
    error('Must invoke main.m from current or parent directory')
end
addpath(fullfile(basedir,'scripts'));

[~,DPRResults]=DPRSim(basedir,datadir);
[~,RPRResults]=RPRSim(basedir,datadir);
disp('========================================================================')
disp('----------------------------------DPR Results---------------------------')
disp('========================================================================')
disp(DPRResults)
disp('')
disp('========================================================================')
disp('----------------------------------RPR Results---------------------------')
disp('========================================================================')
disp(RPRResults)