function [vh,vinhomo,y1]=reveal_real_measurements(N,path,datapath,els,skipcurr,skipvolt,name,perturbations)
%homogeneous data
filename1=[datapath '\data_mat_files\datamat_1_0'];
load(filename1);
fprintf('Select a phantom test\n')
sel1=perturbations(1);
while ~ismember(sel1,1:8)
    msg='Phantom type (1-8):';
    sel1=input(msg);
end
possibles=zeros(8,7);
possibles(1,1:4)=1:4;
possibles(2,1:6)=1:6;
possibles(3,1:6)=1:6;
possibles(4,1:4)=1:4;
possibles(5,1:2)=1:2;
possibles(6,1:7)=1:7;
possibles(7,1:2)=1:2;
possibles(8,1:6)=1:6;
fprintf('Select a phantom object\n')
sel2=perturbations(2);
while ~ismember(sel2,possibles(sel1,:))||sel2==0
    msg='Phantom objects :';
    sel2=input(msg);
end
switch skipcurr
    case 0 %adjacent
        lim1=1;
        lim2=16;
    case 1
        lim1=17;
        lim2=32;
    case 2
        lim1=33;
        lim2=48;
    case 3
        lim1=49;
        lim2=64;
    otherwise
        error('An unexpected error occurred\n')
end
vhm=Uel(:,lim1:lim2);
vhm=reshape(vhm,256,1);
vhm=vhm.*els;
vhm=vhm(vhm~=0);
vh.meas=vhm;
clear CurrentPattern MeasPattern Uel
filename2=[datapath '\data_mat_files\datamat_',num2str(sel1),'_',num2str(sel2)];
load(filename2);
vinhomo=Uel(:,lim1:lim2);
vinhomo=reshape(vinhomo,256,1);
vinhomo=vinhomo.*els;
vinhomo=vinhomo(vinhomo~=0);
filename3=[datapath '\target_photos\target_photos\fantom_',num2str(sel1),'_',num2str(sel2),'.jpg'];
if exist(filename3)~=0
    [y1,y2]=imread(filename3,'jpg');
else
    y1='image';
end

end