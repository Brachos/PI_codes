% function [] = climb(MTOW,WING)

% Conversion    
    feet = 3.28084;  % m to ft
    lbf  = 0.224809; % N to lbf
    kt   = 1.943844; % m/s to kt
    
% Wing geometry
    Sw = WING.Sw;
    AR = WING.A;
    CD_0 = 0.019; % Obtained with drag study
    e = 1.78*(1-0.045*AR^0.68) - 0.64;
    k = 1/(e*pi*AR);

% Available thrust
    % h = (36306 : 1 : 36310)/feet; % [m] Altitude -> h_max = 36307 ft
    % h = (35000 : 1 : 40000)/feet; % [m] Altitude 
    % h = (0 : 10000 : 30000)/feet; % [m] Altitude 
    h = [0 5000 10000 20000 30000 350000]/feet; % [m] Altitude
    
    h1 = 'At sea level';
    h2 = 'At 5,000 ft';
    h3 = 'At 10,000 ft';
    h4 = 'At 20,000 ft';
    h5 = 'At 30,000 ft';
    h6 = 'At 35,000 ft';
    
    [rho, a, T, P] = atmos(h);
    [~,~,~,P0] = atmos(0);

    A = zeros(1,length(h));
    Z = zeros(1,length(h));
    X = zeros(1,length(h));

% Static conditions
    T0 = 13580; % [N] Static thrust
    BPR = 4.1;
    G = 1.1;
    g = 9.80665; % [m/s^2]
    W = MTOW*g;
    
    for i = 1 : length(h)
        A(i) = -0.4327*(P(i)/P0)^2 + 1.3855*(P(i)/P0) + 0.0427;
        X(i) = 0.1377*(P(i)/P0)^2 - 0.4374*(P(i)/P0) + 1.3003;
        Z(i) = 0.9106*(P(i)/P0)^2 - 1.7736*(P(i)/P0) + 1.8697;
    end

    M = 0 : 0.001 : 2;
    Thrust = zeros(length(h),length(M));

    for i = 1 : length(h)
        for j = 1 : length(M)
            Thrust(i,j) = T0*(A(i) - 0.377*(1+BPR)*Z(i)/sqrt(G*(1+0.82*BPR))*P(i)/P0*M(j) + (0.23+0.19*sqrt(BPR))*X(i)*P(i)/P0*M(j)^2);
        end
    end
    
    gamma = 3*pi/180; % [deg] climb angle
    
%     gamma = (0:3:180)*pi/180;
%     h_max = zeros(1,length(gamma));
%     
%     for zaza = 1 : length(gamma)
    
    q  = zeros(length(h),length(M));
    CL = zeros(length(h),length(M));
    CD = zeros(length(h),length(M));
    D  = zeros(length(h),length(M));
    Vv = zeros(length(h),length(M));
    
    for i = 1 : length(h)
        for j = 1 : length(M)
            V(i,j)  = M(j)*a(i);
            q(i,j)  = 0.5*rho(i)*(V(i,j))^2;
            CL(i,j) = W*cos(gamma)/(q(i,j)*Sw);
            CD(i,j) = CD_0 + k*CL(i,j)^2;
            D(i,j)  = q(i,j)*Sw*CD(i,j);
            Vv(i,j) = (Thrust(i,j)-D(i,j))*V(i,j)/W;
        end
    end
    
    figure
    hold on
    for i = 1 : length(h)
        plot(V(i,:)*kt,Vv(i,:)*feet*60)
%         plot(V(i,:)/a(i),Vv(i,:)*feet*60)
        [~,max_pos] = max(Vv(i,:));
        p2 = plot(V(i,max_pos)*kt,Vv(i,max_pos)*feet*60,'x','Color',[0 112 127]/256)
    end
    
    hold off
    legend(h1,h2,h3,h4,h5,h6)
    legend(p2,'Optimal point')
    xlabel('True airspeed [KTSA]')
    ylabel('Rate of climb [fpm]')
    ylim([0 max(Vv(1,:))*60*feet])
      
% % Maximum possible altitude
%     tol = 0.1;
%     for i = 1 : length(h)
%         Max = max(Vv(i,:));
%         if abs(Max)< tol
%             num = i;
%             break
%         end
%     end
%     
%     h_max(zaza) = h(num)*feet;
%     end
%     
%     figure
%     plot(gamma*180/pi,h_max)
%     xlabel('climb angle [deg]')
%     ylabel('Altitude [ft]')