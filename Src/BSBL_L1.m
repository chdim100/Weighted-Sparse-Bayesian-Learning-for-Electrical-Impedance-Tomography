function [sigmaSBL1]=BSBL_L1(xc,yc,Jpixel,dv,h,thetamax,bhta)
%
fprintf('L1- Sparse Bayesian Learning (SBL-L1):\n')
%
M=size(Jpixel,1); L=size(Jpixel,2);
%%%%% no of clusters
g=L-h+1;
Correlations=retrieve_pattern_coupling(xc,yc,h);
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
zi=ones(g,1);
thresholdgamma=0.1;
%%%weights(?)
%wi=0.001*rand(g,1);
go=0.01*sum(sqrt(1/(M-1)*sum(abs(dv-mean(dv,1)).^2,1)),2);
%
So=zeros(g*h,g*h);
Sinv=zeros(g*h,g*h);
%
for i=1:g
    Bi{i}=toeplitz(0.9.^(0:h-1));
    So(h*(i-1)+1:i*h,1:g*h)=sparse([zeros(size(Bi{i},1),(i-1)*h) gi(i)*Bi{i} zeros(size(Bi{i},1),h*(g-1)-(i-1)*h)]);
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
    Sx=So-SoPhiT*(Sv\Phi*So);
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
            gi(i)=1/sqrt(zi(i))*sqrt(muxiT*Biblock\muxi);
        else
            gi(i)=giold(i);
        end
        zi=trace(Biblock*Phi_iT*Sv\Phi_i);
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
    end
    
    %%%%update go using (23)
    go=1/M*(norm(dv-Phi*mu)^2+summary);
    Bi=Binew;
    epsilon=norm(mu-muxprev)/norm(mu);
    fprintf('Iteration %2.0f gives error of %2.4f\n',theta,epsilon)
    theta=theta+1;
end
sigmaSBL1=Psi*mu;

% figure
% plot(dv,'b','LineWidth',2)
% hold on
% plot(Jpixel*sigmaSBL1,'r','LineWidth',2)
% legend({'\delta v','J\cdot\delta\sigma'})
% title('\delta v vs J\cdot\delta\sigma')
%
figure
scatter3(xc,yc,sigmaSBL1,125,sigmaSBL1,'filled')
view([0 90])
title('B-SBL')
colormap jet
colorbar
end