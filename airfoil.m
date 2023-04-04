c_root = 1.4888;
c_tip  = 0.3*c_root;
span   = 6.7740;
sweep  = 15;

coord  = dlmread('SC(2)-0714.txt');
c_upper = coord(1:103,:);
c_lower = coord(end:-1:104,:);
cont = [c_lower;c_upper];

p_root = c_root*coord;
p_tip  = c_tip*coord;
R_root = [cosd(2.5) -sind(2.5);sind(2.5) cosd(2.5)]; % Rotation matrix
R_tip  = [cosd(0.5) -sind(0.5);sind(0.5) cosd(0.5)]; % Rotation matrix

% Reduction pf the number of points
p_root(10:2:100,:) = [];
p_root(end-10:-2:end-80,:)= [];
p_root(8:2:52,:) = [];
p_root(44:2:end-8,:) = [];
p_tip(10:2:100,:) = [];
p_tip(end-10:-2:end-80,:)= [];
p_tip(8:2:52,:) = [];
p_tip(44:2:end-8,:) = [];

% Add AOA
p_root = p_root*R_root;
p_tip  = p_tip*R_tip;

% Add sweep angle
dx_tip = c_root/4 + tand(sweep)*span/2 - c_tip/4;
p_tip(:,1) = p_tip(:,1) + dx_tip;

% Add y coordinate
p_tip_left(:,1)  = p_tip(:,1);
p_tip_right(:,1) = p_tip(:,1);
p_tip_left(:,2)  = -0.5*span*ones(1,length(p_tip(:,1)));
p_tip_right(:,2) = 0.5*span*ones(1,length(p_tip(:,1)));
p_tip_left(:,3)  = p_tip(:,2);
p_tip_right(:,3) = p_tip(:,2);

p_root(:,3) = zeros(1,length(p_tip(:,1)));
p_root = p_root(:,[1 3 2]);

% Plot
c_upper    = p_root(1:34,:);
c_lower    = p_root(35:end,:);
curve_root = [c_upper;c_lower(end:-1:1,:)];
c_upper    = p_tip_left(1:34,:);
c_lower    = p_tip_left(35:end,:);
curve_tip_left  = [c_upper;c_lower(end:-1:1,:)];
c_upper    = p_tip_right(1:34,:);
c_lower    = p_tip_right(35:end,:);
curve_tip_right = [c_upper;c_lower(end:-1:1,:)];

travers = zeros(length(curve_root(:,1)),3,3);

for i = 1 : length(curve_root(:,1))
    travers(i,:,:) = [curve_tip_left(i,:);curve_root(i,:);curve_tip_right(i,:)];
end
% leading  = [curve_tip_left(1,:);curve_root(1,:);curve_tip_right(1,:)];
% trailing = [curve_tip_left(34,:);curve_root(34,:);curve_tip_right(34,:)];

figure
plot3(curve_root(:,1),curve_root(:,2),curve_root(:,3))
hold on
plot3(curve_tip_left(:,1),curve_tip_left(:,2),curve_tip_left(:,3))
plot3(curve_tip_right(:,1),curve_tip_right(:,2),curve_tip_right(:,3))
for i = 1 : length(curve_root(:,1))
    plot3(travers(i,:,1),travers(i,:,2),travers(i,:,3))
end
% plot3(leading(:,1),leading(:,2),leading(:,3))
% plot3(trailing(:,1),trailing(:,2),trailing(:,3))
axis equal
xlim([-1 2]);
ylim([-4 4]);
zlim([-0.2 0.2]);

% writematrix(curve_root,'coord_root.txt');
% writematrix(curve_tip_left,'coord_tip_left.txt');
% writematrix(curve_tip_right,'coord_tip_right.txt');