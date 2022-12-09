%%%%%%%%%%%%%%%%%% Source Code (Main) for Weighted Bound-Optimization Block
%%%%%%%%%%%%%%%%%% Sparse Bayesian Learning EIT. Comparisons with Single-Step
%%%%%%%%%%%%%%%%%% and Iterative Tikhonov Regularization, Total Variation
%%%%%%%%%%%%%%%%%% Regularization and Simple EM-Based BSBL. Quantitative
%%%%%%%%%%%%%%%%%% metrics. 

%%%%%%Cite as: 
%%%%%% Dimas, Christos, Vassilis Alimisis, and Paul P. Sotiriadis. 
%%%%%%"Electrical Impedance Tomography using a Weighted Bound-Optimization 
%%%%%% Block Sparse Bayesian Learning Approach" 
%%%%% 2022 IEEE 22th International Conference on Bioinformatics and Bioengineering (BIBE).
%%%% IEEE, 2022.

%%%%%%%path to experimental data
datapath='Real_Data';
%%%%%%%path to EIDORS
Eidorspath= %%%%%%%%% SET A PATH TO EIDORS !
if ~exist('startuped')||startuped==0
    try 
        run ([Eidorspath,'eidors-v3.9-ng\eidors\startup.m'])
        startuped=1;
    catch 
        error('This code needs the EIDORS package to run properly, please download from https://eidors3d.sourceforge.net/download.shtml')
    end
end
%%%%% select current skip pattern skipcurr and voltage skip pattern
%%%%% skipvolt
skipcurr=0;
skipvolt=0;
%%%%% number of electrodes
N=16;
clc
close all

%%%%%%examples ---general forms
%%%%% A) perturabances=[x0,y0,R,shape,conductivity_change]
%%%%% x0,y0: center of perturbation
%%%%% R: radius of perturbation
%%%%% shape: possible choices: 1 for cycle, 2 for semicycle y>0,
%%%%% 3 for semicycle x>0, 4 for semicycle x<0, 5 for semicycle y<0,
%%%%% 6 for annulus
%%%%% n lines for the above form for n perturbations
%%%%% B) perturbances=[n m] for experimental data from the University of
%%%%% Eastern Finland. Source: V. K. A. Hauptmann and S. Siltanen,
%%%%% “Open 2D Electrical Impedance Tomography data archive,” 2017,
%%%%% arXiv:1704.01178.
%%%%% each pair n,m represents the corresponding fantom (see
%%%%% Real_Data/targets.jpg)
%%%%% C) if perturbances is scalar, then in-vivo data from EIDORS and
%%%%% V. Guardo et al., “A superheterodyne serial data acquisition system for-
%%%%% Electrical impedance tomography,” in Proc. 15th Annu. Int. Conf. IEEE
%%%%% Eng. Med. Biol. Soc., 1993, pp. 86–87.  is loaded, consisting a number
%%%%% of EIT thoracic imaging frames for a single breath cycle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% examples for perturbances

 perturbances=[0.4 -0.2 0.3 6 -0.35;...
     -0.4 -0.2 0.3 6 -0.35;...
     0 0.4 0.3 6 -0.35];

% perturbances=[0.4 0 0.2 1 -0.45;...
%     -0.4 0 0.2 1 -0.45;...
%     0 -0.4 0.3 5 0.35];
% load('reference_1.mat')

% perturbances=[0.4 0 0.2 1 -0.3;...
%     -0.4 0 0.2 1 -0.3;...
%     0 -0.4 0.2 1 0.3;...
%     0 0.4 0.2 1 0.3];

%  perturbances=[0.4 0.4 0.2 1 -0.3;...
%      0.4 -0.4 0.2 1 0.3;...
%      -0.4 -0.4 0.2 1 -0.3;...
%      -0.4 0.4 0.2 1 0.3];
%  load('reference_2.mat')

%perturbances=[3 4];

%perturbances=[2 1];

%%%%check sim or real
if size(perturbances,2)>2 %%%case simulated measurements
    %%%%%% Add noise to the simulated measurements, by defining an SNR
    SNR=40;
    sim=1;
    %%%%%% injected current amplitude
    I=1e-03;
    sim_model_code='j2d3c';
    %%%create simulation model and get the measurements
    [sim_model2,vh,vi,dv]=Create_2D_Simulation_Circular(skipcurr,skipvolt,...
        N,I,SNR,sim_model_code,perturbances);
elseif size(perturbances,2)==2 %%% case real measurements
    I=2e-03;
    [stim,els] = mk_stim_patterns(N,1,[0,skipcurr+1],[0,skipvolt+1],{},I);
    [vh,vi,y1]=reveal_real_measurements(N,'path',datapath,els,skipcurr,skipvolt,'xaxa',perturbances);
    dv=-(vi-vh.meas);
    sim=0;
    clear sim_model sim_model2
elseif size(perturbances,2)==1 %%%case thoracic measurements
    I=2e-03;
    [stim,els] = mk_stim_patterns(N,1,[0,skipcurr+1],[0,skipvolt+1],{},I);
    load montreal_data_1995
    dvt=double(zc_resp(els,15)-zc_resp(els,1));
    vh.meas=double(zc_resp(els,1));
    vh.meas=vh.meas/max(abs(vh.meas));
    dv=dvt; dv=dv/1000; clear sim_model sim_model2
else
    error('Inappropriate size of matrix "perturbances"');
end

%%%set pixelwise model resolution factor
if size(perturbances,2)~=1
    rec_model_code='d2c3';
    coarseres=41;
else
    rec_model_code='c2t3';
    coarseres=0.3;
end
%%%create (pixelwise) reconstruction model
%
inv_model=mk_common_model(rec_model_code,N);
[coarse_fwd_mdl]=make_2d_geom(inv_model,N,coarseres);
[imagemodel_pixelized,xc,yc,non_uniform_elems,elementsorted1D,elementsorted2D]...
    =make_inverse_geometry(N,skipcurr,skipvolt,coarse_fwd_mdl,I);
Vho_str=fwd_solve(imagemodel_pixelized); Vho=Vho_str.meas;
%%%%estimate background (~=1 in simulated cases)
sref=1/abs((Vho'*vh.meas)/(Vho'*Vho));
%
LFEM=length(imagemodel_pixelized.elem_data);
uniform_elems=setdiff(1:LFEM,non_uniform_elems);
% Calculate Jacobian matrix (pixelwise domain)
[Jfull,Jpixel]=...
    calculate_pure_Jacobian(imagemodel_pixelized,...
    non_uniform_elems,elementsorted1D,sref);
%
L=size(Jpixel,2);

%%%%% get the ground truth
if ~exist('sim_model2') %%%case real
    [ImCond,xflatten,yflatten]=color_to_conductivity_decoder(perturbances,datapath);
    path='Quick_Data';
    reference=Map_image_to_pixel([xflatten yflatten],ImCond,xc,yc,path);
else %%%case simulated measurements
    reference=get_ground_truth(sim_model2,sim_model_code,xc,yc);
end
%%%%%
%Tikhonov hyperparameters
if sim==1
    lambdaT=1e-05;
    lambdaN=25*1e-01;
    lambdaL=2*1e-06;
elseif size(perturbances,2)==2 %%%case real experiments
    lambdaT=3*1e-01;
    lambdaN=25*1e-01;
    lambdaL=5*1e-01;
else  %%%%in-vivo
    lambdaT=20;
    lambdaN=50;
    lambdaL=10;
    sim=0;
end

%%%%%%Single-Step Tikhonov
dsigmaT=single_step_GTR(Jpixel,xc,yc,dv,'Standard Tikhonov',lambdaT,sim);
%%%%%%Single-Step NOSER
dsigmaN=single_step_GTR(Jpixel,xc,yc,dv,'NOSER',lambdaN,sim);
%%%%%%Single-Step Laplace
dsigmaL=single_step_GTR(Jpixel,xc,yc,dv,'Laplace',lambdaL,sim);


%%%%%%Iterative Approaches
%Iterative Tikhonov
elim=0.5*1e-10;
dsigmaITT=iterative_GTR(imagemodel_pixelized,xc,yc,Jpixel,dv,lambdaT,...
    'Standard Tikhonov',non_uniform_elems,uniform_elems,elementsorted1D,...
    elementsorted2D,elim,sref,sim);
%Iterative NOSER
dsigmaITN=iterative_GTR(imagemodel_pixelized,xc,yc,Jpixel,dv,lambdaN,...
    'NOSER',non_uniform_elems,uniform_elems,elementsorted1D,elementsorted2D,elim,sref,sim);
%Iterative Laplace
dsigmaITL=iterative_GTR(imagemodel_pixelized,xc,yc,Jpixel,dv,lambdaL,...
    'Laplace',non_uniform_elems,uniform_elems,elementsorted1D,elementsorted2D,elim,sref,sim);
%%%%%%%
%%%%%%%

%%%%%%Iterative TV
%
TV=1;
if TV==1
    cofactor=1.2;
    if sim==1
        lambda=0.6*1e-06;
        bhta=cofactor*1e-03;
    elseif size(perturbances,2)==2 %%%case real experiments
        lambda=1e-02;
        bhta=cofactor*1e-06;
    else
        lambda=1e-02;
        bhta=cofactor*1e-06;
    end
    elim=0.5*1e-10; maxIT=20;
    dsigmaTV=Total_Var(imagemodel_pixelized,xc,yc,Jpixel,dv,lambda,bhta,cofactor,...
        non_uniform_elems,uniform_elems,elementsorted1D,elementsorted2D,elim,maxIT,sref,sim);
    %%%%
end
% select if Sparse Bayesian Learning will be performed
SBL=1;
if SBL==1
    %%%%%%%%%%%%%%%%%%%%%%%%SPARSE BAYESIAN LEARNING
    %%%%%cluster size
    h=4;
    %%%% pattern coupling parameter
    bhta=1;
    %%%% maximum number of iterations
    thetamax=6;
    std=max(diff(xc));
    %%%%Weighted-Convex (MM) Non-Recursive
    tic
    sigmaSBLw=BSBL_weighted(xc,yc,Jpixel,dv,h,thetamax,dsigmaL,std,sim);
    toc
    %%%%%MM-Recursive
    if sim==0
        dvp=dv/max(abs(dv));
    else
        dvp=dv;
    end
    tic
    sigmaSBL=Block_Sparse_Bayesian_Learning(xc,yc,Jpixel,dvp,h,thetamax-1,bhta,sim);
    toc
    %  
    %%%%%%% metrics: CC Correlation Coefficient and RRE Relative
    %%%%%%% Reconstruction Error
    rr1=corrcoef(sigmaSBLw,reference);
    ccWBOBSBL=rr1(1,2);
    norm_reference=(reference-1)/max(abs(reference-1));
    norm_sigmaSBLw=sigmaSBLw/max(abs(sigmaSBLw));
    RRE_WBOBSBL=norm(norm_sigmaSBLw-norm_reference)/norm(norm_reference);
    
    %
    rr2=corrcoef(sigmaSBL,reference);
    ccBSBL=rr2(1,2);
    norm_sigmaSBL=sigmaSBL/max(abs(sigmaSBL));
    RRE_BSBL=norm(norm_sigmaSBL-norm_reference)/norm(norm_reference);
    
    %
    rr3=corrcoef(dsigmaITL,reference);
    ccIGTR=rr3(1,2);
    norm_dsigmaITL=dsigmaITL/max(abs(dsigmaITL));
    RRE_IGTR=norm(dsigmaITL-norm_reference)/norm(norm_reference);
    
    rr4=corrcoef(dsigmaTV,reference);
    ccTV=rr4(1,2);
    norm_dsigmaTV=dsigmaTV/max(abs(dsigmaTV));
    
    RRE_TV=norm(norm_dsigmaTV-norm_reference)/norm(norm_reference);
    
end
%
