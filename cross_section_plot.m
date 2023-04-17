%% 2D plot 

%Fuselage skin 

figure1 = figure(1) ;clf;set(figure1,'defaulttextinterpreter','latex');
hold on
plot(Fuselage.y(1, :), Fuselage.z(1, :),'color', '#00707F', 'LineWidth', 2); 
plot(Fuselage.y(1, :), -Fuselage.z(1, :),'color', '#00707F', 'LineWidth', 2);
plot(0, 0,'.', 'Color', 'k','MarkerSize',40);                                   %center
axis equal
%Stringers
plot(Stringer.position_y(1, :), Stringer.position_z(1,:),'.', 'Color','#00707F','MarkerSize',40)
xlabel("$y$ [m] ")
ylabel("$z$ [m] ")
box on
set(gca,'FontSize',20)
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
set(gca,'XMinorTick','off','YMinorTick','off')
hgexport(figure1,'fuselage_cross_section_2D');
%% 3D plot 

figure2 = figure ;
set(figure2,'defaulttextinterpreter','latex');
hold on
%Stringers
plot3(Stringer.position_x(1, :), Stringer.position_y(1, :), Stringer.position_z(1, :) ,'.', 'Color', '#00707F','MarkerSize',20)
plot3(Stringer.position_x(2, :), Stringer.position_y(2, :), Stringer.position_z(2, :) ,'.', 'Color', '#00707F','MarkerSize',20)
plot3(Stringer.position_x(3, :), Stringer.position_y(3, :), Stringer.position_z(3, :) ,'.', 'Color', '#00707F','MarkerSize',20)

 
for i=1:Stringer.nb
    plot3([Stringer.position_x(1, i) Stringer.position_x(3, i)], [Stringer.position_y(1, i) Stringer.position_y(3, i)], [Stringer.position_z(1, i) Stringer.position_z(3, i)],'color','#00707F', 'LineWidth', 2 )   
end 
%cross-section
plot3(Fuselage.x(1, :), Fuselage.y(1, :),Fuselage.z(1, :),'color', '#00707F', 'LineWidth', 2); 
plot3(Fuselage.x(1, :), Fuselage.y(1, :),-Fuselage.z(1, :),'color', '#00707F', 'LineWidth', 2);
plot3(Fuselage.x(2, :), Fuselage.y(2, :),Fuselage.z(2, :),'color', '#00707F', 'LineWidth', 2); 
plot3(Fuselage.x(2, :), Fuselage.y(2, :),-Fuselage.z(2, :),'color', '#00707F', 'LineWidth', 2);
plot3(Fuselage.x(3, :), Fuselage.y(3, :),Fuselage.z(3, :),'color', '#00707F', 'LineWidth', 2); 
plot3(Fuselage.x(3, :), Fuselage.y(3, :),-Fuselage.z(3, :),'color', '#00707F', 'LineWidth', 2);

xlabel("$x$ [m] ")
ylabel("$y$ [m] ")
zlabel("$z$ [m] ")
box on
view([-20 20])
set(gca,'FontSize',20)
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
set(gca,'XMinorTick','off','YMinorTick','off', 'ZMinorTick','off')
hgexport(figure2,'fuselage_cross_section_3D');