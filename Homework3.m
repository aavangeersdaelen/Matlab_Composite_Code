% This script calculates the ABD Matrix and midplane and ply strains from given
% lamina properties; calculates Tsai-Wu Failure Index for Each Lamina
% Author: Colin Bumgarner
% Last Updated 4/2/19
% Update: Corrected stress/strain calculations; added multiplier to moment
% to induce failure
clear
close all
format short

% Define Material Properties
% In-Plane Elastic Constants
E1 = 30e6;  %Psi
E2 = .75e6;   %Psi
NU12 = 0.25;
G12 = .375e6;  %Psi

% Hygrothermal Properties
alpha_1 = -1.7e-7; %1/*C
alpha_2 = 1.6e-5; %1/*C
beta_1 = 0;
beta_2 = 0;

% Strength Properties
Xt = 150e3; %Psi
Yt = 6e3; %Psi
S = 10e3; %Psi
Xc = -100e3; %Psi
Yc = -17e3; %Psi

% Define Ply
theta = [0,90,45,-45,-45,45,90,0];
% If laminate is symmetric use true
sym = true;
num_lamina = length(theta);
t_ply = .25/length(theta);    %in
t = [1,1,1,1,1,1,1,1]*t_ply; %ply thickness array
h_laminate = num_lamina*t_ply;

Delta = 1 - NU12^2*E2/E1;
Q = [[E1 NU12*E2 0 ]; [NU12*E2 E2 0 ]; [0 0 G12*Delta]]/Delta;

% Define Loading Conditions
% If given Force & Moment Resultants use Type 1
% If given midplane strain & curvatures use Type 2
Type = 1;
% Force & Moment Resultants
Nx = 1000; %lbf/in
Ny = 1000; %lbf/in
Nxy = 1000; %lbf/in

Mx = -100; %lbf
My = -100; %lbf
Mxy = -100; %lbf

% Moment Scaling Factor q
q = 1.92;

% Midplane Strains & Curvatures
epsilon_x = 0; 
epsilon_y = 0;
epsilon_xy = 0;

kappa_x = 0;
kappa_y = 0;
kappa_xy = 0;

% Hygrothermal Effects
delta_t = 0; %*F
delta_c = 0; 

% Initialize Matrices before loop
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
N_thermal = zeros(3,1);
N_moisture = zeros(3,1);

zb = -h_laminate/2; % Bottom of lamina k=1
if Type == 1
    for k = 1:num_lamina
        % Note: For simplicity of using transpose relations, this method
        % actually uses T-inverse as T
        m = cosd(-theta(k));
        n = sind(-theta(k));
        T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, (m^2 - n^2)];
    
        Qbar = T*Q*(T');
    
        zk_mp = zb + t(k)/2; % mid plane of lamina k
        zb = zb + t(k); %bottom of next lamina
        A = A + Qbar*t(k);
        if sym == true
            B = zeros(3,3);
        else
            B = B - Qbar*t(k)*zk_mp;
        end
        D = D + Qbar*(t(k)*zk_mp^2+(t(k)^3)/12);
        
        % Environmental Loading Resultants
        alpha_global = T'^-1*[alpha_1;alpha_2;0];
        N_thermal = N_thermal + (Qbar*(alpha_global)*(delta_t)*(t(k)));
        beta_global = T'^-1*[beta_1;beta_2;0];
        N_moisture = N_moisture + (Qbar*(beta_global)*(delta_c)*(t(k)));
        
    end
  
    % Calculate Mid-Plane Total Strains
    ABD = [[A B];[B D]];

    % Mechanical Force & Moment Resultants
    N = [Nx;Ny;Nxy];
    M = [Mx;My;Mxy];
    
    % Moment scaling factor q
    M = q*M;
    
    % Total Force & Moment Resultants
    N = N + N_thermal + N_moisture;
    
    
    eps0 = A\N;     % Mid-Plane Strains
    kappa = D\M;    % Mid-Plane Curvatures

    % Calculation of individual ply strains & Stresses (global)
    % Need to calculate strains for top and bottom of each lamina
    % re-initialize zb
    zb = -h_laminate/2;

    % Initialize matrices for individual ply data storage
    temp1 = length(theta);
    ply_bottom_strains = zeros(3,1,num_lamina);
    ply_middle_strains = zeros(3,1,num_lamina);
    ply_top_strains = zeros(3,1,num_lamina);
    ply_bottom_stresses = zeros(3,1,num_lamina);
    ply_middle_stresses = zeros(3,1,num_lamina);
    ply_top_stresses = zeros(3,1,num_lamina);
    ply_bottom_strains_local = zeros(3,1,num_lamina);
    ply_middle_strains_local = zeros(3,1,num_lamina);
    ply_top_strains_local = zeros(3,1,num_lamina);
    ply_bottom_stresses_local = zeros(3,1,num_lamina);
    ply_middle_stresses_local = zeros(3,1,num_lamina);
    ply_top_stresses_local = zeros(3,1,num_lamina);
    tsai_wu_index_bottom = zeros(1,num_lamina);
    tsai_wu_index_middle = zeros(1,num_lamina);
    tsai_wu_index_top = zeros(1,num_lamina);
    mlaminate = zeros(3,2*num_lamina+1);
    zb_stored = zeros(1,num_lamina);
    zk_mp_stored = zeros(1,num_lamina);
    zt_stored = zeros(1,num_lamina);

    for x = 1:temp1(1)
        m = cosd(theta(x));
        n = sind(theta(x));
        T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, (m^2 - n^2)];
    
        Qbar = (T^-1)*Q*(T'^-1);
        
        zt = zb + t(x); % top of lamina k
        zk_mp = zb + t(k)/2; % mid plane of lamina k
        
        zb_stored(:,x)= zb;
        zk_mp_stored(:,x) = zk_mp;
        zt_stored(:,x) = zt;
        
        alpha_global = T^-1*[alpha_1;alpha_2;0];
        beta_global = T^-1*[beta_1;beta_2;0];
        
        epsb = eps0 + zb*kappa - alpha_global*delta_t - beta_global*delta_c;
        epsm = eps0 + zk_mp*kappa - alpha_global*delta_t - beta_global*delta_c;
        epst = eps0 + zt*kappa - alpha_global*delta_t - beta_global*delta_c;
        ply_bottom_strains(:,:,x) = epsb;
        ply_middle_strains(:,:,x) = epsm;
        ply_top_strains(:,:,x) = epst;
        zb = zb + t(x); %bottom of next lamina
        
        ply_bottom_stresses(:,:,x) = Qbar*(ply_bottom_strains(:,:,x));
        ply_middle_stresses(:,:,x) = Qbar*(ply_middle_strains(:,:,x));
        ply_top_stresses(:,:,x) = Qbar*(ply_top_strains(:,:,x));

    end
else
    %Finish code for Type 2: Given midplane strains & Curvatures Calculate
    %Force & Moment Resultants
    
end

% Calculation of Tsai-Wu Failure Index
% Strength Tensor Terms
F1 = (1/Xt)+(1/Xc);
F11 = -1/(Xt*Xc);
F2 = (1/Yt)+(1/Yc);
F22 = -1/(Yt*Yc);
F6 = 0;
F66 = (1/S^2);
F12 = 0; %Biaxial Stress Assumption

% Calculation of Individual Ply Failure Idices
for x = 1:temp1(1)
     m = cosd(theta(x));
     n = sind(theta(x));
     T = [m^2, n^2, 2*m*n; n^2, m^2, -2*m*n; -m*n, m*n, (m^2 - n^2)];
    
     Qbar = (T^-1)*Q*(T'^-1);
     
     R=[1 0 0;0 1 0;0 0 2];
     
     %Transform to local lamina c.s.
     ply_bottom_stresses_local(:,:,x) = T*ply_bottom_stresses(:,:,x);
     ply_middle_stresses_local(:,:,x) = T*ply_middle_stresses(:,:,x);
     ply_top_stresses_local(:,:,x) = T*ply_top_stresses(:,:,x);
     ply_bottom_strains_local(:,:,x) = R*T*(R^-1)*ply_bottom_strains(:,:,x);
     ply_middle_strains_local(:,:,x) = R*T*(R^-1)*ply_middle_strains(:,:,x);
     ply_top_strains_local(:,:,x) = R*T*(R^-1)*ply_top_strains(:,:,x);

     
     % Temporary Store of Local Stresses for easier calculation typing
     % bottom
     sigma_1b = ply_bottom_stresses_local(1,:,x);
     sigma_2b = ply_bottom_stresses_local(2,:,x);
     sigma_12b = ply_bottom_stresses_local(3,:,x);
     % middle
     sigma_1m = ply_middle_stresses_local(1,:,x);
     sigma_2m = ply_middle_stresses_local(2,:,x);
     sigma_12m = ply_middle_stresses_local(3,:,x);
     % top
     sigma_1t = ply_top_stresses_local(1,:,x);
     sigma_2t = ply_top_stresses_local(2,:,x);
     sigma_12t = ply_top_stresses_local(3,:,x);
     
     % Calculate Tsai-Wu Failure Index
     % bottom
     tsai_wu_index_bottom(:,x) = F1*sigma_1b+F2*sigma_2b+F6*sigma_12b+F11*sigma_1b^2+...
         F22*sigma_2b^2+F66*sigma_12b^2+2*F12*sigma_1b*sigma_2b;
     % middle
     tsai_wu_index_middle(:,x) = F1*sigma_1m+F2*sigma_2m+F6*sigma_12m+F11*sigma_1m^2+...
         F22*sigma_2m^2+F66*sigma_12m^2+2*F12*sigma_1m*sigma_2m;
     % top
     tsai_wu_index_top(:,x) = F1*sigma_1t+F2*sigma_2t+F6*sigma_12t+F11*sigma_1t^2+...
         F22*sigma_2t^2+F66*sigma_12t^2+2*F12*sigma_1t*sigma_2t;
end

% % Display Results
% fprintf('--------------------------------------------------------------\n')
% fprintf('Material Properties\n');
% fprintf('\n')
% fprintf('Elastic Constants & Hygrothermal Properties \n');
% fprintf('E1(PSI)= %4.0f \n',E1);
% fprintf('E2(PSI)= %4.0f \n',E2);
% fprintf('NU12 = %4.2f \n',NU12);
% fprintf('G12(PSI)= %4.0f \n',G12);
% fprintf('Alpha 1 = %7e \n',alpha_1);
% fprintf('Alpha 2 = %7e \n',alpha_2);
% fprintf('\n')
% fprintf('Strength Properties \n');
% fprintf('Xt(PSI)= %4.0f \n',Xt);
% fprintf('Yt(PSI)= %4.0f \n',Yt);
% fprintf('Xc(PSI)= %4.0f \n',Xc);
% fprintf('Yc(PSI)= %4.0f \n',Yc);
% fprintf('S(PSI)= %4.0f \n',S);
% fprintf('\n')
% fprintf('Applied Force & Moment Resultants \n');
% fprintf('Applied Nx (lbf/in)= %4.2f \n',N(1));
% fprintf('Applied Ny (lbf/in)= %4.2f \n',N(2));
% fprintf('Applied Nxy (lbf/in)= %4.2f \n',N(3));
% fprintf('Applied Mx (lbf/in)= %4.2f \n',Mx);
% fprintf('Applied My (lbf/in)= %4.2f \n',My);
% fprintf('Applied Mxy (lbf/in)= %4.2f \n',Mxy);
% fprintf('\n')
% fprintf('Environmental Effects \n');
% fprintf('Temperature Delta = %7e \n',delta_t);
% fprintf('--------------------------------------------------------------\n')
% 
% A
% B
% D
% fprintf('--------------------------------------------------------------\n')
% 
% for k = 1:num_lamina
%     ply_count = ['Bottom of Ply No. ', num2str(k)];
%     disp(ply_count)
%     fprintf('Z = %4.6f \n', zb_stored(:,k))
%     fprintf('Ply Angle (Degrees)= %4.2f\n',theta(k))
%     fprintf('Epsilon 1 = %8.6f \n',ply_bottom_strains_local(1,:,k));
%     fprintf('Epsilon 2 = %8.6f \n',ply_bottom_strains_local(2,:,k));
%     fprintf('Epsilon 12 = %8.6f \n',ply_bottom_strains_local(3,:,k));
%     fprintf('Epsilon X = %8.6f \n',ply_bottom_strains(1,:,k));
%     fprintf('Epsilon Y = %8.6f \n',ply_bottom_strains(2,:,k));
%     fprintf('Epsilon XY = %8.6f \n',ply_bottom_strains(3,:,k));
%     fprintf('\n')
%     fprintf('Sigma X (PSI)= %4.2f \n',ply_bottom_stresses(1,:,k));
%     fprintf('Sigma Y (PSI)= %4.2f \n',ply_bottom_stresses(2,:,k));
%     fprintf('Sigma XY (PSI)= %4.2f \n',ply_bottom_stresses(3,:,k));
%     fprintf('Sigma 1 (PSI)= %4.2f \n',ply_bottom_stresses_local(1,:,k));
%     fprintf('Sigma 2 (PSI)= %4.2f \n',ply_bottom_stresses_local(2,:,k));
%     fprintf('Sigma 12 (PSI)= %4.2f \n',ply_bottom_stresses_local(3,:,k));
%     fprintf('Tsai-Wu Failure Index = %4.6f \n',tsai_wu_index_bottom(k));
%     fprintf('\n')
%     fprintf('--------------------------------------------------------------\n')
%     fprintf('--------------------------------------------------------------\n')
%     fprintf('\n')
%     ply_count = ['Middle of Ply No. ', num2str(k)];
%     disp(ply_count)
%     fprintf('Z = %4.6f \n', zk_mp_stored(:,k))
%     fprintf('Ply Angle (Degrees)= %4.2f\n',theta(k))
%     fprintf('Epsilon 1 = %8.6f \n',ply_middle_strains_local(1,:,k));
%     fprintf('Epsilon 2 = %8.6f \n',ply_middle_strains_local(2,:,k));
%     fprintf('Epsilon 12 = %8.6f \n',ply_middle_strains_local(3,:,k));
%     fprintf('Epsilon X = %8.6f \n',ply_middle_strains(1,:,k));
%     fprintf('Epsilon Y = %8.6f \n',ply_middle_strains(2,:,k));
%     fprintf('Epsilon XY = %8.6f \n',ply_middle_strains(3,:,k));
%     fprintf('\n')
%     fprintf('Sigma X (PSI)= %4.2f \n',ply_middle_stresses(1,:,k));
%     fprintf('Sigma Y (PSI)= %4.2f \n',ply_middle_stresses(2,:,k));
%     fprintf('Sigma XY (PSI)= %4.2f \n',ply_middle_stresses(3,:,k));
%     fprintf('Sigma 1 (PSI)= %4.2f \n',ply_middle_stresses_local(1,:,k));
%     fprintf('Sigma 2 (PSI)= %4.2f \n',ply_middle_stresses_local(2,:,k));
%     fprintf('Sigma 12 (PSI)= %4.2f \n',ply_middle_stresses_local(3,:,k));
%     fprintf('Tsai-Wu Failure Index = %4.6f \n',tsai_wu_index_middle(k));
%     fprintf('\n')
%     fprintf('--------------------------------------------------------------\n')
%     fprintf('\n')
%     ply_count = ['Top of Ply No. ', num2str(k)];
%     disp(ply_count)
%     fprintf('Z = %4.6f \n', zt_stored(:,k))
%     fprintf('Ply Angle (Degrees)= %4.2f \n',theta(k))
%     fprintf('Epsilon 1 = %8.6f \n',ply_top_strains_local(1,:,k));
%     fprintf('Epsilon 2 = %8.6f \n',ply_top_strains_local(2,:,k));
%     fprintf('Epsilon 12 = %8.6f \n',ply_top_strains_local(3,:,k));
%     fprintf('Epsilon X = %8.6f \n',ply_top_strains(1,:,k));
%     fprintf('Epsilon Y = %8.6f \n',ply_top_strains(2,:,k));
%     fprintf('Epsilon XY = %8.6f \n',ply_top_strains(3,:,k));
%     fprintf('\n')
%     fprintf('Sigma X (PSI)= %4.2f \n',ply_top_stresses(1,:,k));
%     fprintf('Sigma Y (PSI)= %4.2f \n',ply_top_stresses(2,:,k));
%     fprintf('Sigma XY (PSI)= %4.2f \n',ply_top_stresses(3,:,k));
%     fprintf('Sigma 1 (PSI)= %4.2f \n',ply_top_stresses_local(1,:,k));
%     fprintf('Sigma 2 (PSI)= %4.2f \n',ply_top_stresses_local(2,:,k));
%     fprintf('Sigma 12 (PSI)= %4.2f \n',ply_top_stresses_local(3,:,k));
%     fprintf('Tsai-Wu Failure Index = %4.6f \n',tsai_wu_index_top(k));
%     fprintf('\n')
%     fprintf('--------------------------------------------------------------\n')
%     fprintf('\n')
%     
%     
% 
% end    
% 
% fprintf('q multiplier value = %4.2f \n',q);