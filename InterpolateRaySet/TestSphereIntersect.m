center = [0 0 0];
radius = 1;
r0 = [0, 0, -2];
k = [0 0 1];
[p,lambda] = SphereIntersect(center, radius, r0, k)
%%
center = [0 0 0];
radius = 1;
r0 = [0, 0, 200];
k = [0 0 1];
[p,lambda] = SphereIntersect(center, radius, r0, k)
%%

center = [0 0 0];
radius = 1;
r0 = [0, 0.99, -2];
k = [0 0 1];
[p,lambda] = SphereIntersect(center, radius, r0, k)
%%
center = [1 2 3];
radius = 2;
r0 = [2,1,400];
k = center + radius * [0 1 0] - r0;
[p,lambda] = SphereIntersect(center, radius, r0, k)


%%

%% plot a sphere
rad = 2;

theta = 0:pi/50:pi;
st = sin(theta);
ct = cos(theta);

phi = 0 : pi/6 : 2 * pi;
figure(1);
clf;
hold on;
for iphi = phi
    xx = rad * st * cos(iphi);
    yy = rad * st * sin(iphi);
    zz = rad * ct;
    plot3(xx,yy,zz,'k');
end

phi = 0:pi/100:2*pi;
for itheta = 0:pi/6:pi
    r = rad * sin(itheta);
    plot3(sin(phi)*r,cos(phi)*r,rad * cos(itheta)*ones(size(phi)),'k');
end

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
% plot a point and its projection

pt = [1,2,-0.5];
%pt = [0,0,3];
po = [0,0,-rad];
loc = rad * pt(1:2) / (rad+pt(3));
scatter3(pt(1),pt(2),pt(3));
scatter3(loc(1), loc(2), 0);
plot3([po(1),pt(1)],[po(2),pt(2)],[po(3),pt(3)]);
plot3([po(1),loc(1)],[po(2),loc(2)],[po(3),0]);

ps = SphereIntersect([0,0,0], rad, po, pt-po);
scatter3(ps(1),ps(2),ps(3));
plot3([po(1),ps(1)],[po(2),ps(2)],[po(3),ps(3)]);

[ex,ey] = RotDir(loc,radius);
len = 0.3* rad;
plot3([ps(1),ps(1)+ex(1)*len],[ps(2),ps(2)+ex(2)*len],[ps(3),ps(3)+ex(3)*len],'g');
plot3([ps(1),ps(1)+ey(1)*len],[ps(2),ps(2)+ey(2)*len],[ps(3),ps(3)+ey(3)*len],'r');
