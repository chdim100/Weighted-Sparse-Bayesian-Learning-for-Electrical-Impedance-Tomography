function [sim_model2,vh,vi,dv]=Create_2D_Simulation_Circular(skipcurr,skipvolt,N,I,SNR,model_code,perturbances)
%homogeneous model
pathway=['Models\' model_code '\'];
%
if ~exist([pathway model_code '.mat'])
    sim_model=mk_image(mk_common_model(model_code,N),1);
else
    load ([pathway model_code '.mat']);
end
[stim,els] = mk_stim_patterns(N,1,[0,skipcurr+1],[0,skipvolt+1],{},I);
sim_model.fwd_model.stimulation = stim;
sim_model.fwd_solve.get_all_meas = 1;
%add perturbation to the model
sim_model2=sim_model;
for object=1:size(perturbances,1)
    X=perturbances(object,1);
    Y=perturbances(object,2);
    R=perturbances(object,3);
    Sc=perturbances(object,4);
    C=perturbances(object,5);
    %X=0.5; Y=0.2; R=0.15; C=-0.05;
    syms x y z
    if Sc==1
        fun=@(x,y,z)(x-X).^2+(y-Y).^2+0*z.^2<R^2;
        sim_model2.elem_data = sim_model2.elem_data + ...
            C*elem_select(sim_model2.fwd_model, fun);
    elseif Sc==2
        %%%% y>0 semicycle
        ctrs= interp_mesh(sim_model2.fwd_model);
        xe= ctrs(:,1); ye= ctrs(:,2);
        re= sqrt((xe-X).^2+(ye-Y).^2);
        block=(ye-Y>0 & re<R);
        sim_model2.elem_data(block)=1+C;
    elseif Sc==3
        %%%% x>0 semicycle
        ctrs= interp_mesh(sim_model2.fwd_model);
        xe= ctrs(:,1); ye= ctrs(:,2);
        re= sqrt((xe-X).^2+(ye-Y).^2);
        block=(xe-X>0 & re<R);
        sim_model2.elem_data(block)=1+C;
    elseif Sc==4
        %%%% x<0 semicycle
        ctrs= interp_mesh(sim_model2.fwd_model);
        xe= ctrs(:,1); ye= ctrs(:,2);
        re= sqrt((xe-X).^2+(ye-Y).^2);
        block=(xe-X<0 & re<R);
        sim_model2.elem_data(block)=1+C;
    elseif Sc==5
        %%%% y<0 semicycle
        ctrs= interp_mesh(sim_model2.fwd_model);
        xe= ctrs(:,1); ye= ctrs(:,2);
        re= sqrt((xe-X).^2+(ye-Y).^2);
        block=(ye-Y<0 & re<R);
        sim_model2.elem_data(block)=1+C;
    elseif Sc==6
        %%%% annulus
        ctrs= interp_mesh(sim_model2.fwd_model);
        xe= ctrs(:,1); ye= ctrs(:,2);
        re= sqrt((xe-X).^2+(ye-Y).^2);
        block=(re<R & re>=R*0.6);
        sim_model2.elem_data(block)=1+C;
    else
        error('Non defined perturbation shape (1-6)')
    end
end
%%%figure the model
figure
down=-0.4; up=0.4; step=0.2;
sim_model3=sim_model2;
sim_model3.elem_data=sim_model3.elem_data-1;
H1=show_fem(sim_model3,1);
colormap jet
% colorb = colorbar;
% caxis([down up])
% colorb.Ticks=down:step:up;
% colorb.TickLabels=down:step:up;
% colorb.TickLabelInterpreter='latex';
% colorb.Label.String = '$\delta\sigma$ $(S/m)$';
% colorb.FontSize=14;
% colorb.Label.FontSize=16;
% colorb.Label.Interpreter = 'latex';
axis off
set(H1, 'edgecolor', [0.5 0.5 0.5]);
set(H1, 'edgealpha', 0.25);

%%%simulate the measurements
vh=fwd_solve(sim_model);
vi=fwd_solve(sim_model2);
%%%add noise to each measurement vector and get the difference
vhn=awgn(vh.meas,SNR,'measured');
vin=awgn(vi.meas,SNR,'measured');
dv=vin-vhn;
%
end