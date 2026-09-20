function I = getImages(i)
    %Outputs the ith image
    global RunProps CineFile CropWindow ROI
    idx = RunProps.firstImage+i-1;
    I = double(ReadCineFileImage(CineFile,idx,0));

    if sum(CropWindow) ~= 0
        I=imcrop(I,CropWindow);
    end

    if sum(ROI(:)) ~= 0
        I(~ROI)=0;
    end
end