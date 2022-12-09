function [sigmaSBLw]=BSBL_weighted(xc,yc,Jpixel,dv,h,thetamax,ds,std,sim)
%
fprintf('weighted Sparse Bayesian Learning (SBL-weighted):\n')
%
M=size(Jpixel,1); L=size(Jpixel,2);
%%%%% no of clusters
g=L-h+1;
%Correlations=retrieve_pattern_coupling(xc,yc,h);
%
Psi=zeros(L,h*g);
%%%%% Inputs
for i=1:g
    Psi(:,(i-1)*h+1:i*h)=sparse([zeros(i-1,h)' eye(h,h)' zeros(L-i-h+1,h)']');
end
%
Phi=Jpixel*Psi;
%
%
%%%%%Initialize
theta=1;
mu=zeros(g*h,1);
gi=ones(g,1);
thresholdgamma=5*1e-02;
%%%weights
%w=rand(g,1);
%w=[zeros(floor(g/2),1)+1; ones(g-floor(g/2),1)];
w=estimate_weights(xc,yc,g,h,ds,std);
go=0.01*sum(sqrt(1/(M-1)*sum(abs(dv-mean(dv,1)).^2,1)),2);
%
So=zeros(g*h,g*h);
Sinv=zeros(g*h,g*h);
%
for i=1:g
    Bi{i}=toeplitz(0.9.^(0:h-1));
    So(h*(i-1)+1:i*h,1:g*h)=sparse([zeros(size(Bi{i},1),(i-1)*h) gi(i)*Bi{i} zeros(size(Bi{i},1),h*(g-1)-(i-1)*h)]);
    Sinv(h*(i-1)+1:i*h,1:g*h)=sparse([zeros(size(Bi{i},1),(i-1)*h) inv(gi(i)*Bi{i}) zeros(size(Bi{i},1),h*(g-1)-(i-1)*h)]);
end
Bitilde=Bi;
Im=eye(M,M);
epsilon=1;
emin=1e-02;
%
while epsilon>emin && theta<=thetamax
    SoPhiT=So*Phi';
    PhiSoPhiT=Phi*SoPhiT;
    %%%%update mux using (15)
    Sv=single(go*Im+PhiSoPhiT);
    muxprev=mu;
    mu=SoPhiT/Sv*dv;
    %%%%update Sx using (16)
    Sx=So-SoPhiT*(Sv\Phi*So);
    %%%%update gi using (22)
    summary=0;
    giold=gi;
    for i=1:g
        Phi_i=single(Phi(:,(i-1)*h+1:i*h));
        Phi_iT=Phi_i';
        Sxi=Sx((i-1)*h+1:i*h,(i-1)*h+1:i*h);
        Biblock=single(squeeze(Bi{i}));
        muxi=mu((i-1)*h+1:i*h);
        muxiT=muxi';
        %%%%%%update g with pattern coupling
        if giold(i)>thresholdgamma
            Ai=muxiT*(Biblock\muxi);
            Bhi=trace(Sv\Phi_i*Biblock*Phi_iT);
            ai=1/w(i);
            bi=w(i);
            zi=(2*ai+sqrt(4*ai^2+4*Bhi*(Ai+2*bi)))/(2*(Ai+2*bi));
            gi(i)=1/zi;
        else
            gi(i)=giold(i);
        end
        summary=summary+...
            trace(Sxi*Phi_iT*Phi_i);
        %%%%update Bi using (18)-(21)
        %%%%update Bitilde (18)
        Bitildei=squeeze(Bitilde{i});
        Bitilde{i}=Bitildei+1/gi(i)*...
            (Sxi+muxi*muxiT);
        Bitildei=squeeze(Bitilde{i});
        %%%compute rtildei (21)
        rtildei=mean(diag(Bitildei,1))/mean(diag(Bitildei));
        %%%compute ri (20)
        ri=sign(rtildei)*min([norm(rtildei) 0.99]);
        Biblocknew=toeplitz(ri.^(0:1:h-1));
        Binew{i}=Biblocknew;
        %%%update So
        So(h*(i-1)+1:i*h,1:g*h)=sparse([zeros(size(Biblocknew,1),(i-1)*h) gi(i)*Biblocknew zeros(size(Biblocknew,1),h*(g-1)-(i-1)*h)]);
        Sinv(h*(i-1)+1:i*h,1:g*h)=sparse([zeros(size(Biblocknew,1),(i-1)*h) inv(gi(i)*Biblocknew) zeros(size(Biblocknew,1),h*(g-1)-(i-1)*h)]);
        
    end
    %%%%update go using (23)
    go=1/M*(norm(dv-Phi*mu)^2+summary);
    Bi=Binew;
    epsilon=norm(mu-muxprev)/norm(mu);
    fprintf('Iteration %2.0f gives error of %2.4f\n',theta,epsilon)
    theta=theta+1;
end
sigmaSBLw=Psi*mu;

% %%%%%%
% figure
% plot(dv,'b','LineWidth',2)
% hold on
% plot(Jpixel*sigmaSBLw,'r','LineWidth',2)
% legend({'\delta v','J\cdot\delta\sigma'})
% title('\delta v vs J\cdot\delta\sigma')
% %%%%%%
%
if sim==1
    down=-0.4; up=0.4; step=0.2;
else
    down=-2e-03; up=2e-03; step=1e-03;
end
pp=get_recplot(xc,yc,sigmaSBLw,down,up,step);
title('Weighted-BSBL')

end