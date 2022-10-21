function [cl,cp] = panel_method(NACA,AOA,V_inf,c)

    N_pannels = 1000;

    epsilon = NACA(1)/100;
    p       = NACA(2)/10;
    tau     = (NACA(3)*10+NACA(4))/100;

    Xi = linspace(0,2*pi,N_pannels+1);
    x  = c/2 * (cos(Xi)+1);


    %% Thickness

    T = zeros(size(x));
    for i = 1:N_pannels
        T(i) = 10*tau*c*(0.2969*sqrt(x(i)/c) - 0.126*(x(i)/c) - 0.3537*(x(i)/c)^2 + 0.2843*(x(i)/c)^3 - 0.1015*(x(i)/c)^4);
    end


    %% Camber

    Y = zeros(size(x));
    for i = 1:N_pannels
        if (x(i)/c) >= 0 && (x(i)/c) <= p
            Y(i) = epsilon*x(i)/p^2 * (2*p - x(i)/c);
        elseif (x(i)/c) >= p && (x(i)/c) <= 1
            Y(i) = epsilon*(c-x(i))/(1-p)^2 * (1 + x(i)/c - 2*p);
        else
            error('x/c fail (camber)');
        end
    end


    %% Theta for coordinates

    theta_coord = zeros(size(x));
    for i = 1:N_pannels
        if (x(i)/c) >= 0 && (x(i)/c) <= p
            theta_coord(i) = atan(2*epsilon/p - 2*epsilon*x(i)/(p^2*c));
        elseif (x(i)/c) >= p && (x(i)/c) <= 1
            theta_coord(i) = atan(-epsilon/(1-p)^2 + epsilon/(1-p)^2 - 2*epsilon*x(i)/(c*(1-p)^2) + 2*epsilon*p/(1-p)^2);
        else
            error('x/c fail (theta)');
        end
    end


    %% Coordinates

    pannels_x = zeros(N_pannels+1,1);
    pannels_y = zeros(N_pannels+1,1);
    for i = 1:N_pannels+1
        if i < floor(N_pannels/2)+1
            pannels_x(i) = x(i)+T(i)*sin(theta_coord(i))/2;
            pannels_y(i) = Y(i)-T(i)*cos(theta_coord(i))/2;
        else
            pannels_x(i) = x(i)-T(i)*sin(theta_coord(i))/2;
            pannels_y(i) = Y(i)+T(i)*cos(theta_coord(i))/2;
        end
    end

    pannel_length = zeros(N_pannels,1);
    pannels_x_mean = zeros(N_pannels,1);
    pannels_y_mean = zeros(N_pannels,1);
    for i = 1:N_pannels
        pannels_x_mean(i) = (pannels_x(i)+pannels_x(i+1))/2;
        pannels_y_mean(i) = (pannels_y(i)+pannels_y(i+1))/2;
        pannel_length(i) = sqrt((pannels_x(i)-pannels_x(i+1))^2 + (pannels_y(i)-pannels_y(i+1))^2);
    end


    %% Theta for matrices

    theta_mat = zeros(size(x));
    for i = 1:N_pannels
        theta_mat(i) = atan2((pannels_y(i+1)-pannels_y(i)),(pannels_x(i+1)-pannels_x(i)));
    end


    %% Matrices

    A = zeros(N_pannels,N_pannels);
    B = zeros(N_pannels,N_pannels);
    C = zeros(N_pannels,N_pannels);
    D = zeros(N_pannels,N_pannels);
    E = zeros(N_pannels,N_pannels);
    F = zeros(N_pannels,N_pannels);
    G = zeros(N_pannels,N_pannels);

    for i = 1:N_pannels
        for j = 1:N_pannels
            A(i,j) = -(pannels_x_mean(i)-pannels_x(j))*cos(theta_mat(j)) - (pannels_y_mean(i)-pannels_y(j))*sin(theta_mat(j));
            B(i,j) = (pannels_x_mean(i)-pannels_x(j))^2 + (pannels_y_mean(i)-pannels_y(j))^2;
            C(i,j) = sin(theta_mat(i)-theta_mat(j));
            D(i,j) = cos(theta_mat(i)-theta_mat(j));
            E(i,j) = (pannels_x_mean(i)-pannels_x(j))*sin(theta_mat(j)) - (pannels_y_mean(i)-pannels_y(j))*cos(theta_mat(j));
            F(i,j) = log10((1 + (pannel_length(j)^2+2*A(i,j)*pannel_length(j)) / B(i,j)))/log10(exp(1));
            G(i,j) = atan(E(i,j)*pannel_length(j)/(A(i,j)*pannel_length(j)+B(i,j)));
        end
    end


    %% Linear system

    % M
    M_1 = zeros(N_pannels,N_pannels);
    M_2 = zeros(N_pannels,1);
    M_3 = zeros(1,N_pannels);
    M_4 = 0;

    for i = 1:N_pannels
        for j = 1:N_pannels
            if i == j
                M_1(i,j) = pi;
                M_2(i) = M_2(i) + 0;
            else
                M_1(i,j) = C(i,j)*F(i,j)/2 - D(i,j)*G(i,j);
                M_2(i) = M_2(i) + (D(i,j)*F(i,j)/2 + C(i,j)*G(i,j));
            end
        end
    end

    for j = 1:N_pannels
        M_3(1,j) = D(N_pannels,j)*F(N_pannels,j)/2 + C(N_pannels,j)*G(N_pannels,j) + D(1,j)*F(1,j)/2 + C(1,j)*G(1,j);
        M_4 = M_4 - (C(N_pannels,j)*F(N_pannels,j)/2 - D(N_pannels,j)*G(N_pannels,j) + C(1,j)*F(1,j)/2 - D(1,j)*G(1,j));
    end

    M = 1/(2*pi) * [M_1 M_2; M_3 M_4];

    % N

    N = V_inf*sin(theta_mat(1)-AOA);
    for i = 2:N_pannels
        N = [N; V_inf*sin(theta_mat(i)-AOA)];
    end

    N = [N; V_inf*(cos(theta_mat(1)-AOA) + cos(theta_mat(N_pannels)-AOA))];

    % Resulting x

    x_res = M\N;

    %% Outputs

    q = x_res(1:N_pannels);

    gamma = x_res(N_pannels+1);
    
    % Lift coefficient

    cl = 2*gamma/(V_inf*c) * sum(pannel_length);

    % Pressure coefficient
    
    V_t = zeros(N_pannels,1);
    for i = 1:N_pannels
        st = 0;
        tt = pi;
        for j = 1:N_pannels
            if i~=j
                st = st + q(j)/(2*pi) * (D(i,j)*F(i,j)/2 + C(i,j)*G(i,j));
                tt = tt + gamma/(2*pi) * (C(i,j)*F(i,j)/2 - D(i,j)*G(i,j));
            end
        end
        V_t(i) = V_inf*cos(theta_mat(i)-AOA) - st + tt;
    end
    
    cp = 1 - (V_t./V_inf).^2;
    
    % Axial and normal forces
    
    L = -trapz(cp(1:N_pannels/2+1)*sin(theta_mat(1:N_pannels/2+1))) + trapz(cp(1:N_pannels/2+1)*sin(theta_mat(1:N_pannels/2+1)));
    D = trapz(cp(1:N_pannels/2+1)*cos(theta_mat(1:N_pannels/2+1))) - trapz(cp(1:N_pannels/2+1)*cos(theta_mat(1:N_pannels/2+1)));
    
    
    %% Plots

    % airfoil
    hold on
    plot(pannels_x,pannels_y)
    scatter(pannels_x_mean,pannels_y_mean,'.')
    axis([0 1 -0.5 0.5])
    hold off
end