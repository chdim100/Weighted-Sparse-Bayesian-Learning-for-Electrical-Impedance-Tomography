function MAP=Map_image_to_pixel(ImageCoords,Imagedata,x,y,path)
L=length(Imagedata);
dir2search=[path,'\MAPS\From_Image_to_2Dpixel\'];
file2search=['Image_',num2str(L),'_to_',num2str(length(x)),'pixels.mat'];
Map_exists=exist([dir2search file2search]);
switch Map_exists
    case 2
        load([dir2search file2search])
        CORRs(CORRs==0)=1;
        MAP=Imagedata(CORRs);
    otherwise
        MAP=zeros(length(x),1);
        CORRs=zeros(length(x),1);
        for point=1:length(x)
            elemin=1;
            mindist=1000;
            r=[x(point) y(point)];
            for element=1:L
                D=sum((ImageCoords(element,:) - r) .^ 2);
                if D<mindist
                    mindist=D;
                    elemin=element;
                end
            end
            MAP(point)=Imagedata(elemin);
            CORRs(point)=elemin;
        end
        file2save=['Image_',num2str(L),'_to_',num2str(length(x)),'pixels.mat'];
        save([dir2search file2save],'CORRs')
end
end