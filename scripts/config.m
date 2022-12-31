function datadir=config(basedir,datadir) 
    if nargin<1
        basedir='..';
    end
    if nargin<2
        datadir=fullfile(basedir,'/data');
    end
    javapath=fullfile(basedir,'bin');
    octavepath=fullfile(basedir,'scripts','octave');

    dpath=javaclasspath();
    if sum(strcmpi(dpath,javapath))>0
        javarmpath(javapath);
    end
    javaaddpath(javapath);

    addpath(octavepath);
end