function [] = createCartPendulumAnimation(X,name,fps)
    %Function that creates an animation from the X matrix obtained from
    %CartPendulum_ReferenceGovernor_HorizontalPosition.m

    v = VideoWriter(name,'Motion JPEG AVI');
    v.Quality = 100;
    v.FrameRate = fps;
    L = 7.5;

    close all
    
    open(v);
    
    for k = 1:size(X,2)
        rr = rectangle('position',[X(1,k)-2 -1 4 2]);
        ll = line([X(1,k) X(1,k)+L*sin(X(3,k))],[0 L*cos(X(3,k))],'color','b','LineWidth',5);
        stt = text(6,-7,sprintf('s = %.2f',X(1,k)));
        xlim([-11.5 11.5])
        ylim([-8 8])
        frame = getframe(gcf);
        writeVideo(v,frame);
        set(rr,'Visible','off')
        set(ll,'Visible','off')
        set(stt,'Visible','off')
    end
    close(v);

end