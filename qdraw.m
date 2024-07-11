function qdraw(A,x,y)
imageq(A)
B=A(x(1):x(2),y(1):y(2));
hold on;
imageql(B,[512-(x(2)-x(1))*2 512],[512-(y(2)-y(1))*2 512])
plot([512-(x(2)-x(1))*2 512-(x(2)-x(1))*2],[512-(y(2)-y(1))*2 512],'r','LineWidth',2);
plot([512-(x(2)-x(1))*2 512],[512-(x(2)-x(1))*2 512-(y(2)-y(1))*2],'r','LineWidth',2);
plot([512 512],[512-(x(2)-x(1))*2 512],'r','LineWidth',2);
plot([512-(x(2)-x(1))*2 512],[512 512],'r','LineWidth',2);

plot([x(1) x(1)],[y(1) y(2)],'r','LineWidth',2);
plot([x(2) x(2)],[y(1) y(2)],'r','LineWidth',2);
plot([x(1) x(2)],[y(1) y(1)],'r','LineWidth',2);
plot([x(1) x(2)],[y(2) y(2)],'r','LineWidth',2);
axis off
end
% plot([0 412],[100 512],'r','LineWidth',2);