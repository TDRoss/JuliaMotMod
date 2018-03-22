tlng = 3000;
figure;
for z = 2:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)));
q.AutoScale='off';
axis([-20000 20000 -20000 20000])
hold on
quiver(filhist(z,:,1)+.5*tlng*cos(filhist(z,:,3)),filhist(z,:,2)+.5*tlng*sin(filhist(z,:,3)),10000*(filhist(z,:,1)-filhist(z-1,:,1)),10000*(filhist(z,:,2)-filhist(z-1,:,2)))
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel(z*.1)
pause(0.1)
end

figure;
writerObj = VideoWriter('AsterFailedLink02Sim.avi');
writerObj.FrameRate = 20;
open(writerObj);
for z = 2:2:size(filhist,1)
    q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)));
q.AutoScale='off';
axis([-20000 20000 -20000 20000])
hold on
quiver(filhist(z,:,1)+.5*tlng*cos(filhist(z,:,3)),filhist(z,:,2)+.5*tlng*sin(filhist(z,:,3)),10000*(filhist(z,:,1)-filhist(z-1,:,1)),10000*(filhist(z,:,2)-filhist(z-1,:,2)))
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel(z*.1)
frame = getframe(gcf);
writeVideo(writerObj, frame);
pause(0.1)
end
close(writerObj);


figure;
writerObj = VideoWriter('Aster3umFilLink9s.avi');
writerObj.FrameRate = 20;
open(writerObj);
for z = 2:2:size(filhist,1)
    q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)));
q.AutoScale='off';
axis([-15000 15000 -6000 6000])
hold on
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel(z*.1)
frame = getframe(gcf);
writeVideo(writerObj, frame);
pause(0.1)
end
close(writerObj);



figure;
writerObj = VideoWriter('AsterLinkSimNoNoise.avi');
writerObj.FrameRate = 20;
open(writerObj);
for z = 2:2:150
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'linewidth',2);
q.AutoScale='off';
axis([-6000 8000 -7000 9000])
hold on
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel(z*.1)
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);

figure;
writerObj = VideoWriter('PlusMinusCollapse01.avi');
writerObj.FrameRate = 20;
open(writerObj);
wp = moptype<3;
wn = ~wp;
for z = 2:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'linewidth',1);
q.AutoScale='off';
axis([-4000 4000 -4000 4000])
hold on
plot(mothist(z,wp,1),mothist(z,wp,2),'*k','markers',1)
plot(mothist(z,wn,1),mothist(z,wn,2),'*m','markers',1)
hold off
xlabel(z*.1)
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);

figure;
writerObj = VideoWriter('PlusMinusSwitch02.avi');
writerObj.FrameRate = 20;
open(writerObj);
wp = moptype<3;
wn = ~wp;
for z = 2:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'linewidth',1);
q.AutoScale='off';
axis([-6000 2000 -5000 5000])
hold on
plot(mothist(z,wp,1),mothist(z,wp,2),'*k','markers',1)
plot(mothist(z,wn,1),mothist(z,wn,2),'*m','markers',1)
hold off
xlabel(z*.1)
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);


figure;
wp = moptype<3;
wn = ~wp;
for z = 2:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'linewidth',1);
q.AutoScale='off';
axis([-10000 10000 -10000 10000])
hold on
plot(mothist(z,wp,1),mothist(z,wp,2),'*k','markers',1)
plot(mothist(z,wn,1),mothist(z,wn,2),'*m','markers',1)
hold off
xlabel(z*.1)
pause(0.1)
end

mting = zeros(size(mforcehist,1),1);
for i = 1:size(mting,1)
    w = mforcehist(i,:)>0;
    mting(i) = mean(mforcehist(i,w));
end
% q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'linewidth',1.5);
% q.AutoScale='off';
% axis([-20000 20000 -20000 20000])
% xlabel(strcat(num2str((z-1)*.1),' (s)'))
% frame = getframe(gcf);
% writeVideo(writerObj, frame);
% pause(0.1)
% end
% close(writerObj);

figure;
for z = 1:2:size(filhist,1)
    xforce = forcehist(z,:,1).*cos(filhist(z,:,3)) - forcehist(z,:,2).*sin(filhist(z,:,3));
    yforce = forcehist(z,:,1).*sin(filhist(z,:,3)) + forcehist(z,:,2).*cos(filhist(z,:,3));
    mag = sqrt(xforce.^2+yforce.^2);
    w = mag>0;
    th = atan2(yforce(w),xforce(w));
q = quiver(filhist(z,w,1)+.5*tlng*cos(filhist(z,w,3)),filhist(z,w,2)+.5*tlng*sin(filhist(z,w,3)),1000*mag(w).*cos(th),1000*mag(w).*sin(th));
q.AutoScale='off';
axis([-10000 10000 -10000 10000])
xlabel(z*.1)
pause(0.1)
end




figure;
writerObj = VideoWriter('AsterFormationWithNoise.avi');
writerObj.FrameRate = 10;
open(writerObj);
for z = 2:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2);
q.AutoScale='off';
axis([-20000 20000 -20000 20000])
hold on
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel('nm')
ylabel('nm')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);


figure;
writerObj = VideoWriter('AsterFormationWithNoiseRedEnds.avi');
writerObj.FrameRate = 10;
open(writerObj);
for z = 2:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-20000 20000 -20000 20000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',10)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel('nm')
ylabel('nm')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);


figure;
writerObj = VideoWriter('MotorA1.0pi_MotorB1.0pi_k40.avi');
writerObj.FrameRate = 10;
open(writerObj);
for z = 1:20:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-18000 18000 -18000 18000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',10)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',1)
hold off
xlabel('nm')
ylabel('nm')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);


figure;
tlng = argvec{9};
for z = 1:5:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-40000 40000 -20000 20000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',1)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',.001)
hold off
xlabel(z*.5)
ylabel('nm')
daspect([1 1 1])
pause(0.1)
end

figure; 
w = filhist(:,1,3)~=0;
sw = sum(w);
for z =1:2:sw
    histogram(mothist(z,:,1),'normalization','probability')
    pause(0.1)  
end


tlng = argvec{9};
figure;
writerObj = VideoWriter('DiskExcitationSim.avi');
writerObj.FrameRate = 10;
open(writerObj);
for z = 1:20:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-25000 25000 -25000 25000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',10)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',.1)
hold off
xlabel('nm')
ylabel('nm')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj);


t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);

tlng = argvec{9};
figure;
for z = 1:5:size(filhist,1)
    h=fill(x,y,'y','LineStyle','none');
    hold on
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*k','markers',3)
plot(mothist(z,:,1),mothist(z,:,2),'.r','markers',.3)
set(h,'FaceAlpha',.3)
hold off
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
pause(0.16)
end

n = 20;
x = linspace(-10,10,n);
y = [x.^2;x;2*x];
p = plot(x,y(1,:),'r', 'LineWidth',5);
% modified jet-colormap
modc = parula(20*4);
cd = [uint8(modc(1:n,:)*255) uint8(ones(n,1))].';
hold on
for i=1:3
 p = plot(x,y(i,:),'r', 'LineWidth',5);  
drawnow
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
end
hold off

modc = cbrewer('div','RdBu',(4));




modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
figure;
set(gcf,'Visible','off')
for z = 1:40:size(filhist,1)
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
%     xps = [0 ; tlng];
%     yps = [0 ;  0];
    plot(xps,yps,'LineWidth',1.5)
    hold on
%     for j = 1:size(xps,2)
%          p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
%         drawnow
%     set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
%     end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
set(gcf,'Visible','on')
% pause(0.16)
end



modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
figure;
writerObj = VideoWriter('DiskExcitationSimFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z =1:10:size(filhist,1)
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
%     xps = [0 ; tlng];
%     yps = [0 ;  0];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)


modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
chng = 0;
figure;
writerObj = VideoWriter('MoveRectangleSimFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
w = find(sum(filhist(:,:,1),2)~=0);
for z =1:2:max(w)
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
%     xps = [0 ; tlng];
%     yps = [0 ;  0];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
if z>=50 && chng == 0
    t = (0.0:.001:1)'*2*pi+t(579);
    t = t(1:end-1);
    x = 15000+15000*cos(t);
    y = 15000*sin(t);
    x = [ -40000; -40000;x;-40000];
    y = [7000; -7000;y;7000];
    w = x>-40000 & x<15000 & y > -7000 & y<7000;
    x = x(~w);
    y = y(~w);
        chng =1;
end
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)


figure;
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
tlng = argvec{9};
writerObj = VideoWriter('MoveSlowDiskSimFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:900
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
if z>30
x = 15000+15000*cos(t) - (z*.5 -15)*60;
end
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)



figure;
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
tlng = argvec{9};
writerObj = VideoWriter('MoveSlowDiskSimFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:900
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-40000 40000 -20000 20000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',1)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',.001)
if z>30
x = 15000+15000*cos(t) - (z*.5 -15)*60;
end
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)
% print( 'LocalAsterSimFrame600','-dpng','-r600');


figure;
tlng = argvec{9};
for z = 1:20:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-40000 40000 -20000 20000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',1)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',.001)
hold off
xlabel(z*.5)
ylabel('nm')
daspect([1 1 1])
pause(0.1)
end


figure;
mag = 10000;
for z = 1:2:size(filhist,1)
   q= quiver(filhist(z,:,1)+mag*cos(filhist(z,:,3)),filhist(z,:,2)+mag*sin(filhist(z,:,3)),mag*forcehist(z,:,1),mag*forcehist(z,:,2));
   q.AutoScale='off';
axis([-40000 40000 -20000 20000])
    daspect([1 1 1])
    pause(0.1)
end

figure;
mag = 10;
for z = 2:10:size(filhist,1)
   q= quiver(filhist(z,:,1)+10000*cos(filhist(z,:,3)),filhist(z,:,2)+10000*sin(filhist(z,:,3)),mag*diff(filhist((z-1):z,:,1),1),mag*diff(filhist((z-1):z,:,2),1));
   q.AutoScale='off';
axis([-40000 40000 -20000 20000])
    daspect([1 1 1])
    set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
    pause(0.1)
end

figure;
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
writerObj = VideoWriter('AsterLocalization_FilVelocity.avi');
writerObj.FrameRate = 10;
open(writerObj);
mag = 10;
for z = 2:10:size(filhist,1)
   q= quiver(filhist(z,:,1)+10000*cos(filhist(z,:,3)),filhist(z,:,2)+10000*sin(filhist(z,:,3)),mag*diff(filhist((z-1):z,:,1),1),mag*diff(filhist((z-1):z,:,2),1));
   hold on
   q.AutoScale='off';
axis([-40000 40000 -20000 20000])
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
    daspect([1 1 1])
set(gca,'xtick',[],'ytick',[])
set(gca,'position',[0 0 1 1],'units','normalized')
box off
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)


figure;
t = (0.0:.001:1)'*2*pi;
x = 15000+15000*cos(t);
y = 15000*sin(t);
writerObj = VideoWriter('AsterLocalization_MotVelocity.avi');
writerObj.FrameRate = 10;
open(writerObj);
mag = 2;
for z = 2:10:size(filhist,1)
    wf = find(mothist(z-1,:,3)>0);
    if size(wf)>0
    f = find(mothist(z-1,wf,3)==mothist(z,wf,3));
    w = wf(f);
   q= quiver(mothist(z,w,1),mothist(z,w,2),mag*diff(mothist((z-1):z,w,1),1),mag*diff(mothist((z-1):z,w,2),1));
   hold on
   q.AutoScale='off';
axis([-40000 40000 -20000 20000])
h=fill(x,y,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
    daspect([1 1 1])
set(gca,'xtick',[],'ytick',[])
set(gca,'position',[0 0 1 1],'units','normalized')
box off
frame = getframe(gcf);
writeVideo(writerObj, frame)
    end
end
close(writerObj)



figure;
tlng = argvec{9};
for z = 1:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-40000 40000 -20000 20000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',1)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',.001)
h=fill(xp,yp,'y','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
xlabel(z*.5)
ylabel('nm')
daspect([1 1 1])
pause(0.1)
end

figure;
t = (0.0:.001:1)'*2*pi;
x = 500*cos(t);
y = 500*sin(t);
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
writerObj = VideoWriter('RectangleFormCoMFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:260
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
h=fill([25000 -25000 -25000 25000],[4000 4000 -4000 -4000],'y','LineStyle','none');
h2=fill(x,y,'r','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
set(h2,'FaceAlpha',.3)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)


xp = [-25000 -25000 -5000 -5000 25000 25000 -5000 -5000];
yp = [15000 -15000 -15000 -4000 -4000 4000 4000 15000];
figure;
t = (0.0:.001:1)'*2*pi;
x = 500*cos(t);
y = 500*sin(t);
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
writerObj = VideoWriter('HammerFormFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:260
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
h=fill(xp,yp,'y','LineStyle','none');
% h2=fill(x,y,'r','LineStyle','none');
hold off
set(h,'FaceAlpha',.15)
% set(h2,'FaceAlpha',.3)
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)


figure;
t = (0.0:.001:1)'*2*pi;
x1 = 25000+8000*cos(t);
y1 = 8000*sin(t);
x2 = -25000+8000*cos(t);
y2 = 8000*sin(t);
lx = [25000 -25000 -25000 25000];
ly =[5000 5000 -5000 -5000];
tlng = argvec{9};
for z = 1:2:size(filhist,1)
q = quiver(filhist(z,:,1),filhist(z,:,2),tlng*cos(filhist(z,:,3)),tlng*sin(filhist(z,:,3)),'LineWidth',2,'ShowArrowHead','off');
q.AutoScale='off';
axis([-40000 40000 -20000 20000])
hold on
plot(filhist(z,:,1)+tlng*cos(filhist(z,:,3)),filhist(z,:,2)+tlng*sin(filhist(z,:,3)),'*r','markers',1)
plot(mothist(z,:,1),mothist(z,:,2),'*k','markers',.001)
fill(x1,y1,'y','LineStyle','none','FaceAlpha',.15);
fill(x2,y2,'y','LineStyle','none','FaceAlpha',.15);
if z>120
fill(lx,ly,'y','LineStyle','none','FaceAlpha',.15);
end
hold off
xlabel(z)
ylabel('nm')
daspect([1 1 1])
pause(0.1)
end


figure;
t = (0.0:.001:1)'*2*pi;
x1 = 25000+8000*cos(t);
y1 = 8000*sin(t);
x2 = -25000+8000*cos(t);
y2 = 8000*sin(t);
lx = [25000 -25000 -25000 25000];
ly =[5000 5000 -5000 -5000];
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
writerObj = VideoWriter('TwoAsterLinkFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:300
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
fill(x1,y1,'y','LineStyle','none','FaceAlpha',.15);
fill(x2,y2,'y','LineStyle','none','FaceAlpha',.15);
if z>120
fill(lx,ly,'y','LineStyle','none','FaceAlpha',.15);
end
hold off
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)



figure;
t = (0.0:.001:1)'*2*pi;
x1 = 25000+8000*cos(t);
y1 = 8000*sin(t);
lx = [25000 -25000 -25000 25000];
ly =[5000 5000 -5000 -5000];
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
writerObj = VideoWriter('OneAsterLinkFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:300
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
fill(x1,y1,'y','LineStyle','none','FaceAlpha',.15);
% fill(x2,y2,'y','LineStyle','none','FaceAlpha',.15);
if z>120
fill(lx,ly,'y','LineStyle','none','FaceAlpha',.15);
end
hold off
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)


figure;
t = (0.0:.001:1)'*2*pi;
x1 = 25000+8000*cos(t);
y1 = 8000*sin(t);
x2 = -25000+8000*cos(t);
y2 = 8000*sin(t);
lx = [25000 -25000 -25000 25000];
ly =[5000 5000 -5000 -5000];
modc = parula(3);
cd = [uint8(modc(1:2,:)*255) uint8(ones(2,1))].';
tlng = argvec{9};
writerObj = VideoWriter('TwoAsterLinkFailFancy.avi');
writerObj.FrameRate = 10;
open(writerObj);
set(gcf,'Visible','off')
for z = 1:5:300
    xps = [filhist(z,:,1) ; filhist(z,:,1)+tlng*cos(filhist(z,:,3))];
    yps = [filhist(z,:,2) ; filhist(z,:,2)+tlng*sin(filhist(z,:,3))];
    for j = 1:size(xps,2)
        if j==2
            set(gcf,'Visible','off')
            hold on
        end
         p = plot(xps(:,j),yps(:,j),'r', 'LineWidth',2);  
        drawnow
    set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    end
plot(mothist(z,:,1),mothist(z,:,2),'.','Color',[.6,0,0.05],'markers',4)
axis([-40000 40000 -15000 15000])
set(gca,'xtick',[],'ytick',[])
daspect([1 1 1])
fill(x1,y1,'y','LineStyle','none','FaceAlpha',.15);
fill(x2,y2,'y','LineStyle','none','FaceAlpha',.15);
if z>120
fill(lx,ly,'y','LineStyle','none','FaceAlpha',.15);
end
hold off
set(gca,'visible','off','position',[0 0 1 1],'units','normalized')
box off
% set(gcf,'Visible','on')
frame = getframe(gcf);
writeVideo(writerObj, frame)
end
close(writerObj)