%Author: LWeissinger, 18.07.2025
%%Assumes that the variable points holds the vertices of a bezier-curve and
%%returns the length of this curve. accuracy improves with increased
%%integer input "gen". "fit" assumes linear interpolated additional
%%vertices (fit aditional vertices per point). This results in the
%%resulting beziercurve fitting the original polygon more (Uncomment the plot line and try for several parameters to compare)


function [distances,distance] = curvedistance(app,points,gen,fit)
h=waitbar(0,sprintf('Initializing distance computation'));
distance=0;
d=size(points,2);
distances=zeros(size(points,1),1);
points_start=zeros(size(points,1)+(size(points,1)-1)*(fit-1),d);
for i=1:(size(points,1)-1)
    for j=1:fit
        if j==1
        points_start(fit*(i-1)+1,:)=points(i,:);
        else
            points_start(fit*(i-1)+j,:)=points(i,:)+(j-1)/fit*(points(i+1,:)-points(i,:));
        end
    end
end
points_start(size(points_start,1),:)=points(size(points,1),:);
s=(size(points_start,1)-1)*gen+1;
newpoints=zeros(s,d);
for run=1:s
    if isvalid(h)
        waitbar (run/s, h,sprintf('Distance computation'))
    end
    t=(run-1)/((size(points_start,1)-1)*gen);
    newpoints(run,:)=deCasteljau(t,points_start);
    if run>1
        step=0;
        for dim=1:size(points_start,2)
            step=step+(newpoints(run,dim)-newpoints(run-1,dim))^2;
        end
        distance=distance+sqrt(step);
    end
    if mod((run-1)/(fit*gen)+1,1)==0
        distances((run-1)/(fit*gen)+1)=distance;
    end
end
cla(app.Display3);
app.Display3.Visible='off';
app.Display2.Visible='on';
cla(app.Display2);
daspect(app.Display2,[1 1 1])
view(app.Display2, [0 0 -1]);
minp=min(points_start,[],1);
maxp=max(points_start,[],1);
switch dim
    case 2
        hold on
        plot(app.Display2,points_start(:,1),points_start(:,2),newpoints(:,1),newpoints(:,2))
        plot(app.Display2,newpoints((((1:size(points,1))-1)*gen*fit+1),1),newpoints((((1:size(points,1))-1)*gen*fit+1),2),'*')
        if minp(1)<maxp(1)
            app.Display2.XLim=[minp(1) maxp(1)]
        end
        if minp(2)<maxp(2)
            app.Display2.YLim=[minp(2) maxp(2)]
        end
        hold off
    case 3
        hold on
        plot3(app.Display2,points_start(:,1),points_start(:,2),points_start(:,3),newpoints(:,1),newpoints(:,2),newpoints(:,3))
        plot3(app.Display2,newpoints((((1:size(points,1))-1)*gen*fit+1),1),newpoints((((1:size(points,1))-1)*gen*fit+1),2),newpoints((((1:size(points,1))-1)*gen*fit+1),3),'*')
        if minp(1)<maxp(1)
            app.Display2.XLim=[minp(1) maxp(1)]
        end
        if minp(2)<maxp(2)
            app.Display2.YLim=[minp(2) maxp(2)]
        end
        if minp(3)<maxp(3)
            app.Display2.ZLim=[minp(3) maxp(3)]
        end
        hold off
    otherwise
        fprintf("no plotable dimension")
end
delete(h);
end

function point = deCasteljau(t, pts)
beta = pts;
    n = size(beta,1);
    for j=1:n
        for k=1:(n-j)
            beta(k,:) = beta(k,:) * (1 - t) + beta(k+1,:) * t;
        end
    end
    point=beta(1,:);
end