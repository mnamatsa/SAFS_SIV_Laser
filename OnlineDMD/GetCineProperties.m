function [cineProps] = GetCineProperties(fileName)
    %Returns some properties for the cine file in filename

    [HRES, cineHandle] = PhNewCineFromFile(fileName);
    if (HRES<0)
	    [message] = PhGetErrorMessage( HRES );
        error(['Cine handle creation error: ' message]);
    end


    pFirstIm = libpointer('int32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_FIRSTIMAGENO, pFirstIm);
    cineProps.firstImage = pFirstIm.Value;
    pImCount = libpointer('uint32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_IMAGECOUNT, pImCount);
    cineProps.numberOfImages = double(pImCount.Value);

    pFramerate = libpointer('int32Ptr',0);
    PhGetCineInfo(cineHandle, PhFileConst.GCI_FRAMERATE, pFramerate);
    cineProps.frameRate = double(pFramerate.Value);    

end