%Program creates an animation from the .mat files created at the end of
%each of the main programs. Uncomment the instructions at the end of each
%program to create these files.
%
% NOTE: Set saveAnimation to 1 to save animetion to AVI file.

close all
clear all
clc

%load('Inverted_Pendulum_on_a_Cart_Stabilization_MPC.mat')
%load('Inverted_Pendulum_on_a_Cart_Movement_MPC.mat')
%load('Inverted_Pendulum_on_a_Cart_Movement_LQR_RG.mat')
%load('Inverted_Pendulum_on_a_Cart_Swing_Up_NMPC.mat')
load('Inverted_Pendulum_on_a_Cart_Swing_Up_Switching_MPC.mat')

%name = "Inverted_Pendulum_on_a_Cart_Stabilization_MPC";
%name = "Inverted_Pendulum_on_a_Cart_Movement_MPC";
%name = "Inverted_Pendulum_on_a_Cart_Movement_LQR_RG";
%name = "Inverted_Pendulum_on_a_Cart_Swing_Up_NMPC";
name = "Inverted_Pendulum_on_a_Cart_Swing_Up_Switching_MPC";

saveAnimation = 0.0;  %Set to 1 to save animation on AVI file.
L = 7.5; %Pendulum length
SwingUp = 1.0; %Set to 1 to give enough vertical space to properly see the swing-up animation.
fps = 10; %Frame per second
Tmax = max(tt);
tt_Vid = 0:(1/fps):Tmax;
anim_steps = length(tt_Vid);
pos_Vid = interp1(tt,pos,tt_Vid,'makima');
ang_Vid = interp1(tt,ang,tt_Vid,'makima');

if saveAnimation >= 0.5
    v = VideoWriter(name,'Motion JPEG AVI');
    v.Quality = 100;
    v.FrameRate = fps;
    open(v);
end

figure(1)
set(gcf, 'color', [1 1 1])
set(gcf, 'position',[100,100,1280,720])
ax = gca;
ax.FontSize = 17; 
set(gca,'TickLabelInterpreter','latex')

grid on
box on

yline(0,'k--','LineWidth',2)
     
for kk = 1:anim_steps
    rr = rectangle('position',[pos_Vid(kk)-2 -0.75 4 1.5], 'Curvature',0.2, 'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0 0],'LineWidth',1);
    ll = rectangle('position',[0 -0.5 0.25 L-0.5], 'Curvature',0.2, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0 0 0],'LineWidth',1);
    ll_g = makehgtform('zrotate', -ang_Vid(kk));
    ll_t_vec = [cos(-ang_Vid(kk)) -sin(-ang_Vid(kk)); sin(-ang_Vid(kk)) cos(-ang_Vid(kk))]*[-0.125; 0.25];
    ll_t = makehgtform('translate', [pos_Vid(kk)+ll_t_vec(1) ll_t_vec(2) 0]);
    ll_gt = hgtransform('Matrix',ll_t*ll_g);
    ll.Parent = ll_gt;
    xlim([-11.5+pos_Vid(kk) 11.5+pos_Vid(kk)])
    xlabel('$s$ (m)','Interpreter','latex','FontSize',16)
    xticks(floor((-10:2:10)+pos_Vid(kk)))
    if SwingUp <= 0.5
        ylim([-2 L+0.5]) %For stabilization and cart movement simulations
    else
        ylim([-L-0.5 L+0.5]) %For swing-up simulations
    end
    pause(1/fps)
    if saveAnimation >= 0.5
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    set(rr,'Visible','off')
    set(ll,'Visible','off')
end

if saveAnimation >= 0.5
    close(v);
end