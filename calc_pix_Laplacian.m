function RtR=calc_pix_Laplacian(x,y,z)
Lx=length(x);
step=max(diff(x));
e=step/10;
Reg=zeros(Lx,Lx);
if nargin<3
    for pixel=1:Lx
       horizontal_neighbors{pixel}=find((y==y(pixel)).*(abs(x-x(pixel))<=(step+e)).*(abs(x-x(pixel))>=(step-e)));
       vertical_neighbors{pixel}=find((x==x(pixel)).*(abs(y-y(pixel))<=(step+e)).*(abs(y-y(pixel))>=(step-e)));
       Lh(pixel)=length(horizontal_neighbors{pixel});
       Lv(pixel)=length(vertical_neighbors{pixel});
       Reg(pixel,pixel)=Lh(pixel)+Lv(pixel);
       Reg(pixel,horizontal_neighbors{pixel})=-1;
       Reg(pixel,vertical_neighbors{pixel})=-1;
    end
else
    for pixel=1:Lx
       horizontal_neighbors{pixel}=find((y==y(pixel)).*z==z(pixel).*(abs(x-x(pixel))<=(step+e)).*(abs(x-x(pixel))>=(step-e)));
       vertical_neighbors{pixel}=find((x==x(pixel)).*z==z(pixel).*(abs(y-y(pixel))<=(step+e)).*(abs(y-y(pixel))>=(step-e)));
       depth_neighbors{pixel}=find((x==x(pixel)).*y==y(pixel).*(abs(z-z(pixel))<=(step+e)).*(abs(z-z(pixel))>=(step-e)));
       Lh(pixel)=length(horizontal_neighbors{pixel});
       Lv(pixel)=length(vertical_neighbors{pixel});
       Ld(pixel)=length(depth_neighbors{pixel});
       Reg(pixel,pixel)=Lh(pixel)+Lv(pixel)+Ld(pixel);
       Reg(pixel,horizontal_neighbors{pixel})=-1;
       Reg(pixel,vertical_neighbors{pixel})=-1;
       Reg(pixel,depth_neighbors{pixel})=-1;
    end
end
RtR=Reg;
if sum(sum(RtR))~=0
    error('Check Laplacian Prior')
end
end