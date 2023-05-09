% function [] = climb(MTOW,WING)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
% Conversion    
    feet = 3.28084;  % m to ft
    lbf  = 0.224809; % N to lbf
    kt   = 1.943844; % m/s to kt
    hp   = 0.00134102; % W to HP
    
% Wing geometry
    Sw = WING.Sw;
    AR = WING.A;
    CD_0 = 0.019; % Obtained with drag study
    e = 1.78*(1-0.045*AR^0.68) - 0.64;
    k = 1/(e*pi*AR);

% Available thrust
    % h = (36306 : 1 : 36310)/feet; % [m] Altitude -> h_max = 36307 ft
    h = (35000 : 1 : 40000)/feet; % [m] Altitude 
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
    
%     for i = 1 : length(h)
%         for j = 1 : length(M)
%             V(i,j)  = M(j)*a(i);
%             q(i,j)  = 0.5*rho(i)*(V(i,j))^2;
%             CL(i,j) = W*cos(gamma)/(q(i,j)*Sw);
%             CD(i,j) = CD_0 + k*CL(i,j)^2;
%             D(i,j)  = q(i,j)*Sw*CD(i,j);
%             Vv(i,j) = (Thrust(i,j)-D(i,j))*V(i,j)/W;
%         end
%     end
%     
%     figure
%     hold on
%     for i = 1 : length(h)
%         plot(V(i,:)*kt,Vv(i,:)*feet*60)
% %         plot(V(i,:)/a(i),Vv(i,:)*feet*60)
%     end
%     for i = 1 : length(h)
%         [~,max_pos] = max(Vv(i,:));
%         plot(V(i,max_pos)*kt,Vv(i,max_pos)*feet*60,'x','Color',[0 112 127]/256)
%     end
%     
%     hold off
%     legend(h1,h2,h3,h4,h5,h6,'$V_{v_{max}}$','interpreter','latex')
% %     legend(p2)
%     xlabel('True airspeed [KTSA]')
%     ylabel('Rate of climb [fpm]')
%     xlim([0 900])
%     ylim([0 max(Vv(1,:))*60*feet])
%     
%     Power_av = Thrust(1,:).*V(1,:);
%     Power_re = D(1,:).*V(1,:);
%     DP = Power_av-Power_re;
%     [~,ind]=max(DP);
%     
%     figure 
%     plot(V(1,:)*kt,Power_av*hp,'Color',[0 112 127]/256)
%     hold on
%     plot(V(1,:)*kt,Power_re*hp,'Color',[240 127 60]/256)
%     plot([V(1,ind)*kt V(1,ind)*kt],[Power_av(ind) Power_re(ind)]*hp,'--k')
%     hold off
%     xlabel('True airspeed [KTAS]','interpreter','latex')
%     ylabel('Available and required power [HP]','interpreter','latex')
%     xlim([0 800])
%     ylim([0 0.6e7*hp])
%     legend('$P_{av}$','$P_{re}$','$\Delta P_{max}$')

    
      
% Maximum possible altitude
    tol = 0.1;
    for i = 1 : length(h)
        Max = max(Vv(i,:));
        if abs(Max)< tol
            num = i;
            break
        end
    end
    
    h_max = h(num)*feet
    
    
   