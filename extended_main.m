%% Extended version of the paper:
% A. Shahmansoori, G. E. Garcia, G. Destino, G. Seco-Granados and H. Wymeersch, 
% "Position and Orientation Estimation Through Millimeter-Wave MIMO in 5G Systems," 
% in IEEE Transactions on Wireless Communications, vol. 17, no. 3, pp. 1822-1835, March 2018.
% This version go with a refinement step using Golden-Section search for
% AOA, AOD and Least Square (LS) for channel amplitude and TOA. It is simplified to work with 1 LOS and 1 NLOS (L = 2) and not work in other cases.
% This extended version was developed by Ngoc-Son Duong, Quoc-Tuan Nguyen, Thai-Mai Dinh-Thi (sondn24@vnu.edu.vn).
clear; close all; clc;
%% System Params
L       = 2;                    % number of paths
Rs      = 100e6;                % total BW in Hz
N       = 10;                   % number of subcarriers
Nt      = 32;                   % number of TX antennas
Nr      = Nt;                   % number of RX antennas
Nb      = Nt;                   % number of beams in dictionary
Ns      = 5;                    % number of beams sent
c       = 3e8;                  % speed of light m/s
Ts      = 1/Rs;                 % sampling period in us
posTx   = [0 0]';               % TX is assumed to be in [0, 0]
posRx   = [4 1]';               % RX's(UE) location
SP      = [2 2];                % Scatter's location (with assume SP is always in the 1st quadrant)
alpha   = 0.2;                  % UE's orientation
h       = 10*ones(1,L);         % channel gain
GR = (sqrt(5) - 1)/2;           % Golden ration
search_range = 0.1;             % search range of Golden-Section search
tol = 5*10^(-5);                % stopping criteria used for Golden-Section search
sigma = 0.01;                   % noise
%% Compute Channel Parameters for L paths
TOA = zeros(1, L); AOD = zeros(1, L); AOA = zeros(1, L);
TOA(1) = norm(posRx)/c;                                                     % LOS TOA
AOD(1) = atan2(posRx(2), posRx(1));                                         % LOS AOD
AOA(1) = atan2(posRx(2), posRx(1)) - pi - alpha;                            % LOS AOA
if(AOA(1) < -pi)
    AOA(1) = 2*pi + AOA(1);
end
for p = 1:L-1
    TOA(p+1) = (norm(SP(p,:)) + norm(posRx - SP(p,:)'))/c;                  % NLOS TOA
    AOD(p+1) = atan2(SP(p,2), SP(p,1));                                     % NLOS AOD
    AOA(p+1) = atan2(SP(p,2) - posRx(2), SP(p,1) - posRx(1)) - alpha;       % NLOS AOA
end
%% Create dictionary
Ut = zeros(Nt,Nb); Ur = zeros(Nr,Nb);
aa = -Nb/2:Nb/2-1; aa = 2*aa/Nb;
for m = 1:Nb
    Ut(:,m) = getResponse(Nt,aa(m))*sqrt(Nt);
    Ur(:,m) = getResponse(Nr,aa(m))*sqrt(Nr);
end
%% Generate channel
H = zeros(Nr,Nt,N); A_rx = zeros(Nr, L); A_tx = zeros(Nt, L); Gamma = zeros(L, L);
for n = 1:N
    for p = 1:L
        A_rx(:,p) = getResponse(Nr,sin(AOA(p)));
        A_tx(:,p) = getResponse(Nt,sin(AOD(p)));
        Gamma(p,p,n) = sqrt(Nr*Nt)*h(p)*exp(-1j*2*pi*TOA(p)*(n-1)/(N*Ts));
        H(:,:,n) = H(:,:,n) + A_rx(:,p)*Gamma(p,p,n)*A_tx(:,p)';
    end
end
%% Generate the observation and beamformers
y = zeros(Nr,Ns,N); signal = zeros(Nr,Ns,N); noise = zeros(Nr,Ns,N); F = zeros(Nt,Ns,N);
for k = 1:Ns
    for n = 1:N
        F(:,k,n) = exp(1j*rand(Nt,1)*2*pi);                                 % random beamformers
        signal(:,k,n) = H(:,:,n)*F(:,k,n);                                  % noise-free signal
        noise(:,k,n) = sigma/sqrt(2)*(randn(Nr,1) + 1i*randn(Nr,1));        % noise
        y(:,k,n) = signal(:,k,n) + noise(:,k,n);                            % received signal with noise
    end
end
%% Vectorize and generation of the basis
yb = zeros(Nr*Ns,N);
Omega = zeros(Nr*Ns,Nb*Nb,N);
for n = 1:N
    yb(:,n) = reshape(y(:,:,n), Nr*Ns,1);
    Omega(:,:,n) = kron((Ut'*F(:,:,n)).',Ur);
end
y_vec = reshape(yb, [Nt*Ns*N 1]);
%% DCS-SOMP (coarse estimate)
[indices, h_hat_val] = DCSSOMP(yb, Omega, L);                               % Perform DCS-SOMP
TOA_hat = zeros(1, L); AOD_hat = zeros(1, L); AOA_hat = zeros(1, L);
ch_amp_hat = mean(abs(h_hat_val),2).';                                      % Coarse channel amplitude
index_vec = (1 : Nb*Nb).';                                                  % Vector of indices used for mapping
index_mat = reshape(index_vec, [Nb, Nb]);                                   % Matrix of indices used for mapping
RC = zeros(L, 2);                                                           % Row-Column pairs used for mapping
for p = 1:L
    [RC(p,1), RC(p,2)] = find(index_mat == indices(p));
end
for p = 1:L
    TOA_hat(p) = -mean(diff(phase(h_hat_val(p,:))))*(N*Ts)/(2*pi);
    if (TOA_hat(p) < 0)
        TOA_hat(p) = TOA_hat(p) + N*Ts;                                     % Coarse TOA
    end
end
if(TOA_hat(1) > TOA_hat(2))                                                 % LOS/NLOS detection
    TOA_hat = sort(TOA_hat);
    tmp = RC(1,:);
    RC(1,:) = RC(2,:);
    RC(2,:) = tmp;
end
for p = 1:L
    AOD_hat(p) = asin(aa(RC(p,2)));                                         % Coarse AOD
    AOA_hat(p) = pi*sign(asin(aa(RC(p,1)))) - asin(aa(RC(p,1)));            % Coarse AOA
end
%% Proposed Refinement
Omega_diag = blkdiag(Omega(:,:,1), Omega(:,:,2), Omega(:,:,3), ...
    Omega(:,:,4), Omega(:,:,5), Omega(:,:,6), Omega(:,:,7), ...
    Omega(:,:,8), Omega(:,:,9), Omega(:,:,10));                             % Default sensing matrix
M = 15;                                                                     % Number of loop
% Variables for plot
LOS_AOA = zeros(1, M); LOS_AOD = zeros(1, M); LOS_TOA = zeros(1, M);
NLOS_AOA = zeros(1, M); NLOS_AOD = zeros(1, M); NLOS_TOA = zeros(1, M);
true_AOA0 = zeros(1,M); true_AOD0 = zeros(1,M); true_TOA0 = zeros(1,M);
true_AOA1 = zeros(1,M); true_AOD1 = zeros(1,M); true_TOA1 = zeros(1,M);
tic
for i = 1:M
    %% Refinement for the LOS path
    H10 = zeros(Nr,Nt,N); Hb10 = zeros(Nb,Nb,N);
    for n = 1:N
        for p = 2
            H10(:,:,n) = H10(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(p)))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
        end
        Hb10(:,:,n) = Ur'*H10(:,:,n)*Ut;
    end
    Hb10_vec = reshape(Hb10, [Nb*Nb*N  1]);
    C0_temp = y_vec - Omega_diag*Hb10_vec/(Nb*Nb);
    % Golden-Section search for LOS AOA
    AOA_hat_t = AOA_hat(1) - search_range;
    AOA_hat_s = AOA_hat(1) + search_range;
    delta_AOA = abs(AOA_hat_t - AOA_hat_s);
    while delta_AOA > tol
        d1 = GR*(AOA_hat_s - AOA_hat_t);
        x1 = AOA_hat_t + d1;
        x2 = AOA_hat_s - d1;
        delta_AOA = x1 - x2;
        H0_x1 = zeros(Nr,Nt,N);
        Hb0_x1 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 1
                H0_x1(:,:,n) = H0_x1(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(x1))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
            end
            Hb0_x1(:,:,n) = Ur'*H0_x1(:,:,n)*Ut;
        end
        Hb0_vec_x1 = reshape(Hb0_x1, [Nb*Nb*N  1]);
        H0_x2 = zeros(Nr,Nt,N);
        Hb0_x2 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 1
                H0_x2(:,:,n) = H0_x2(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(x2))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
            end
            Hb0_x2(:,:,n) = Ur'*H0_x2(:,:,n)*Ut;
        end
        Hb0_vec_x2 = reshape(Hb0_x2, [Nb*Nb*N  1]);
        if(norm(C0_temp - Omega_diag*Hb0_vec_x1/(Nb*Nb))^2 > norm(C0_temp - Omega_diag*Hb0_vec_x2/(Nb*Nb))^2)
            AOA_hat_s = x1;
        else
            AOA_hat_t = x2;
        end
    end
    AOA_hat(1) = (x1 + x2)/2;
    LOS_AOA(i) = AOA_hat(1);
    % Golden-Section search for LOS AOD
    AOD_hat_t = AOD_hat(1) - search_range;
    AOD_hat_s = AOD_hat(1) + search_range;
    delta_AOD = abs(AOD_hat_t - AOD_hat_s);
    while delta_AOD > tol
        d2 = GR*(AOD_hat_s - AOD_hat_t);
        x3 = AOD_hat_t + d2;
        x4 = AOD_hat_s - d2;
        delta_AOD = x3 - x4;
        H0_x3 = zeros(Nr,Nt,N);
        Hb0_x3 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 1
                H0_x3(:,:,n) = H0_x3(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(1)))*sqrt(Nt)*getResponse(Nt,sin(x3))';
            end
            Hb0_x3(:,:,n) = Ur'*H0_x3(:,:,n)*Ut;
        end
        Hb0_vec_x3 = reshape(Hb0_x3, [Nb*Nb*N  1]);
        H0_x4 = zeros(Nr,Nt,N);
        Hb0_x4 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 1
                H0_x4(:,:,n) = H0_x4(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(1)))*sqrt(Nt)*getResponse(Nt,sin(x4))';
            end
            Hb0_x4(:,:,n) = Ur'*H0_x4(:,:,n)*Ut;
        end
        Hb0_vec_x4 = reshape(Hb0_x4, [Nb*Nb*N  1]);
        if(norm(C0_temp - Omega_diag*Hb0_vec_x3/(Nb*Nb))^2 > norm(C0_temp - Omega_diag*Hb0_vec_x4/(Nb*Nb))^2)
            AOD_hat_s = x3;
        else
            AOD_hat_t = x4;
        end
    end
    AOD_hat(1) = (x3 + x4)/2;
    LOS_AOD(i) = AOD_hat(1);
    % LS estimate for channel amplitude and TOA of the LOS path
    ur0 = getResponse(Nr,sin(AOA_hat(1)))*sqrt(Nr);                         % Obtain new LOS AOA-related atom in each step
    ut0 = getResponse(Nt,sin(AOD_hat(1)))*sqrt(Nt);                         % Obtain new NLOS AOA-related atom in each step
    H10re = zeros(Nr,Nt,N);
    Hb10re = zeros(1,1,N);
    for n = 1:N
        Omega_re(:,:,n) = kron((ut0'*F(:,:,n)).',ur0);                      % Construct new basis based on 2 atoms related to current estimate of AOA and AOD;
    end
    Are0 = blkdiag(Omega_re(:,:,1), Omega_re(:,:,2), Omega_re(:,:,3), ...
        Omega_re(:,:,4), Omega_re(:,:,5), Omega_re(:,:,6), Omega_re(:,:,7), ...
        Omega_re(:,:,8), Omega_re(:,:,9), Omega_re(:,:,10));                % New basis for all subcarriers
    for n = 1:N
        for p = 2
            H10re(:,:,n) = H10re(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(p)))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
        end
        Hb10re(:,:,n) = ur0'*H10re(:,:,n)*ut0;
    end
    Hb10_vecre = reshape(Hb10re, [1*1*N  1]);
    C0_temp_re = y_vec - Are0*Hb10_vecre/(Nb*Nb);
    h0_re = inv(Are0'*Are0)*Are0'*C0_temp_re;
    ch_amp_hat(1) = mean(abs(h0_re));
    TOA_hat(1) = -mean(diff(phase(h0_re)))*(N*Ts)/(2*pi);
    LOS_TOA(i) = TOA_hat(1);
    %% Refinement for the NLOS path
    H01 = zeros(Nr,Nt,N); Hb01 = zeros(Nb,Nb,N);
    for n = 1:N
        for p = 1
            H01(:,:,n) = H01(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(p)))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
        end
        Hb01(:,:,n) = Ur'*H01(:,:,n)*Ut;
    end
    Hb01_vec = reshape(Hb01, [Nb*Nb*N  1]);
    C1_temp = y_vec - Omega_diag*Hb01_vec/(Nb*Nb);
    % Golden-Section search for NLOS AOA
    AOA_hat_t2 = AOA_hat(2) - search_range;
    AOA_hat_s2 = AOA_hat(2) + search_range;
    delta_AOA2 = abs(AOA_hat_t2 - AOA_hat_s2);
    while delta_AOA2 > tol
        d5 = GR*(AOA_hat_s2 - AOA_hat_t2);
        x9 = AOA_hat_t2 + d5;
        x10 = AOA_hat_s2 - d5;
        delta_AOA2 = x9 - x10;
        H1_x9 = zeros(Nr,Nt,N); Hb1_x9 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 2
                H1_x9(:,:,n) = H1_x9(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(x9))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
            end
            Hb1_x9(:,:,n) = Ur'*H1_x9(:,:,n)*Ut;
        end
        Hb1_vec_x9 = reshape(Hb1_x9, [Nb*Nb*N  1]);
        H1_x10 = zeros(Nr,Nt,N); Hb1_x10 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 2
                H1_x10(:,:,n) = H1_x10(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(x10))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
            end
            Hb1_x10(:,:,n) = Ur'*H1_x10(:,:,n)*Ut;
        end
        Hb1_vec_x10 = reshape(Hb1_x10, [Nb*Nb*N  1]);
        if(norm(C1_temp - Omega_diag*Hb1_vec_x9/(Nb*Nb))^2 > norm(C1_temp - Omega_diag*Hb1_vec_x10/(Nb*Nb))^2)
            AOA_hat_s2 = x9;
        else
            AOA_hat_t2 = x10;
        end
    end
    AOA_hat(2) = (x9 + x10)/2;
    NLOS_AOA(i) = AOA_hat(2);
    % Golden-Section search for NLOS AOD
    AOD_hat_t2 = AOD_hat(2) - search_range;
    AOD_hat_s2 = AOD_hat(2) + search_range;
    delta_AOD2 = abs(AOD_hat_t2 - AOD_hat_s2);
    while delta_AOD2 > tol
        d6 = GR*(AOD_hat_s2 - AOD_hat_t2);
        x11 = AOD_hat_t2 + d6;
        x12 = AOD_hat_s2 - d6;
        delta_AOD2 = x11 - x12;
        H1_x11 = zeros(Nr,Nt,N); Hb1_x11 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 2
                H1_x11(:,:,n) = H1_x11(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(p)))*sqrt(Nt)*getResponse(Nt,sin(x11))';
            end
            Hb1_x11(:,:,n) = Ur'*H1_x11(:,:,n)*Ut;
        end
        Hb1_vec_x11 = reshape(Hb1_x11, [Nb*Nb*N  1]);
        H1_x12 = zeros(Nr,Nt,N); Hb1_x12 = zeros(Nb,Nb,N);
        for n = 1:N
            for p = 2
                H1_x12(:,:,n) = H1_x12(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(p)))*sqrt(Nt)*getResponse(Nt,sin(x12))';
            end
            Hb1_x12(:,:,n) = Ur'*H1_x12(:,:,n)*Ut;
        end
        Hb1_vec_x12 = reshape(Hb1_x12, [Nb*Nb*N  1]);
        if(norm(C1_temp - Omega_diag*Hb1_vec_x11/(Nb*Nb))^2 > norm(C1_temp - Omega_diag*Hb1_vec_x12/(Nb*Nb))^2)
            AOD_hat_s2 = x11;
        else
            AOD_hat_t2 = x12;
        end
    end
    AOD_hat(2) = (x11 + x12)/2;
    NLOS_AOD(i) = AOD_hat(2);
    % LS estimate for channel amplitude and TOA of the NLOS path
    ur1 = getResponse(Nr,sin(AOA_hat(2)))*sqrt(Nr);
    ut1 = getResponse(Nt,sin(AOD_hat(2)))*sqrt(Nt);
    H01re = zeros(Nr,Nt,N);
    Hb01re = zeros(1,1,N);
    for n = 1:N
        Omega_re(:,:,n) = kron((ut1'*F(:,:,n)).',ur1);                      % Construct new basis based on 2 atoms related to current estimate of AOA and AOD;
    end
    Are1 = blkdiag(Omega_re(:,:,1), Omega_re(:,:,2), Omega_re(:,:,3), ...
        Omega_re(:,:,4), Omega_re(:,:,5), Omega_re(:,:,6), Omega_re(:,:,7), ...
        Omega_re(:,:,8), Omega_re(:,:,9), Omega_re(:,:,10));                % New basis for all subcarriers
    for n = 1:N
        for p = 1
            H01re(:,:,n) = H01re(:,:,n) + ch_amp_hat(p)*exp(-1j*2*pi*TOA_hat(p)*(n-1)/(N*Ts))*sqrt(Nr)*getResponse(Nr,sin(AOA_hat(p)))*sqrt(Nt)*getResponse(Nt,sin(AOD_hat(p)))';
        end
        Hb01re(:,:,n) = ur1'*H01re(:,:,n)*ut1;
    end
    Hb01_vecre = reshape(Hb01re, [1*1*N  1]);
    C1_temp_re = y_vec - Are1*Hb01_vecre/(Nb*Nb);
    h1_re = inv(Are1'*Are1)*Are1'*C1_temp_re;
    ch_amp_hat(2) = mean(abs(h1_re));
    TOA_hat(2) = -mean(diff(phase(h1_re)))*(N*Ts)/(2*pi);
    NLOS_TOA(i) = TOA_hat(2);
end
toc
for i = 1:M
    true_AOA0(i) = AOA(1); true_AOD0(i) = AOD(1); true_TOA0(i) = TOA(1); 
    true_AOA1(i) = AOA(2); true_AOD1(i) = AOD(2); true_TOA1(i) = TOA(2);
end
figure(1)
subplot(231); 
plot(LOS_AOA); hold on; grid on; plot(true_AOA0); title('LOS AOA')
subplot(232)
plot(LOS_AOD); hold on; grid on; plot(true_AOD0); title('LOS AOD')
subplot(233)
plot(LOS_TOA); hold on; grid on; plot(true_TOA0); title('LOS TOA')
subplot(234)
plot(NLOS_AOA); hold on; grid on; plot(true_AOA1); title('NLOS AOA')
subplot(235)
plot(NLOS_AOD); hold on; grid on; plot(true_AOD1); title('NLOS AOD')
subplot(236)
plot(NLOS_TOA); hold on; grid on; plot(true_TOA1); title('NLOS TOA')
%% Re-estimating
% Update sensing matrix using refined AOA and AOD;
Ur(:,RC(1,1)) = getResponse(Nr,sin(AOA_hat(1)))*sqrt(Nr);
Ut(:,RC(1,2)) = getResponse(Nt,sin(AOD_hat(1)))*sqrt(Nt);
Ur(:,RC(2,1)) = getResponse(Nr,sin(AOA_hat(2)))*sqrt(Nr);
Ut(:,RC(2,2)) = getResponse(Nt,sin(AOD_hat(2)))*sqrt(Nt);
for n = 1:N
    Omega(:,:,n) = kron((Ut'*F(:,:,n)).',Ur);                               % Refined sensing matrix
end
% Re-estimate channel
[indices_re, h_hat_val_re] = DCSSOMP(yb, Omega, L);
ch_amp_hat = mean(abs(h_hat_val_re),2).';                                   % Refined channel amplitude
for p = 1:L
    TOA_hat(p) = -mean(diff(phase(h_hat_val_re(p,:))))*(N*Ts)/(2*pi);       % Refined time delay (TOA)
end
if(TOA_hat(1) > TOA_hat(2))
    TOA_hat = sort(TOA_hat);
end