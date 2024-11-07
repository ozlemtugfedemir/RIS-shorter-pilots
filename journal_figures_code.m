% Number of BS antennas
Mh = 4;
Mv = 4;
M = Mh*Mv;

% Number of RIS elements
Nh = 16;
Nv = 16;
N = Nh*Nv;

% Antenna (element) spacings in terms of wavelength
dH = [0.25; 0.125];
dV = [0.25; 0.125];


Ncluster = 1; %number of clusters
varphi = zeros(Ncluster,4); % azimuth angles
theta = zeros(Ncluster,4);  % elevation angles

% Nominal angles for R_g' and R_h
varphi(:,1:2) = pi/4;
theta(:,1:2)  = 0;

% Nominal angles for R_g
varphi(:,3) = -pi/4;
theta(:,3) = -pi/6;

% Nominal angles for R_w
varphi(:,4) = pi/5;
theta(:,4) = -pi/12;


ASD_varphi = deg2rad(5); %angular standard deviation for azimuth angle
ASD_theta = deg2rad(5);  %angular standard deviation for elevation angle

MposY = vec(repmat(0:Mh-1,Mv,1)*dH(1));
MposZ = vec(repmat((0:Mv-1).',1,Mh)*dV(1));

NposY = vec(repmat(0:Nh-1,Nv,1)*dH(2));
NposZ = vec(repmat((0:Nv-1).',1,Nh)*dV(2));

RisoBS = zeros(M,M);
RisoRIS = zeros(N,N);

RgP = zeros(M,M);
Rh = zeros(N,N);
Rg = zeros(N,N);
Rw = zeros(N,N);

%Parameters for cosine antenna pattern from the paper [33]
aExp = 1;
bExp0 = 1;
bExp = bExp0+1;

scaleeB = zeros(4,1);
for n = 1:Ncluster
    scaleeB(1) = scaleeB(1) + 1/Ncluster*(cos(varphi(n,1))^aExp*cos(theta(n,1))^bExp);
    scaleeB(2) = scaleeB(2) + 1/Ncluster*(cos(varphi(n,2))^aExp*cos(theta(n,2))^bExp);
    scaleeB(3) = scaleeB(3) + 1/Ncluster*(cos(varphi(n,3))^aExp*cos(theta(n,3))^bExp);
    scaleeB(4) = scaleeB(4) + 1/Ncluster*(cos(varphi(n,4))^aExp*cos(theta(n,4))^bExp);
end
%% Generate spatial correlation matrices for BS
for row = 1:M
    row
    for column = row+1:M

        distanceY = MposY(row)-MposY(column);
        distanceZ = MposZ(row)-MposZ(column);
        RisoBS(row,column) =  sin(2*pi*sqrt(distanceY^2+distanceZ^2))/(2*pi*sqrt(distanceY^2+distanceZ^2));


        for n = 1:Ncluster


            A = exp(1i*2*pi*(distanceY*sin(varphi(n,1))*cos(theta(n,1))+distanceZ*sin(theta(n,1))));
            B = 2*pi*distanceY*cos(varphi(n,1))*cos(theta(n,1));
            C = -2*pi*distanceY*cos(varphi(n,1))*sin(theta(n,1));
            D = -2*pi*distanceY*sin(varphi(n,1))*sin(theta(n,1))+2*pi*distanceZ*cos(theta(n,1));

            sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
            X = cos(varphi(n,1))^aExp*(cos(theta(n,1))^bExp-1i*bExp*cos(theta(n,1))^(bExp-1)*sin(theta(n,1))*ASD_theta^2*D);
            Y = aExp*cos(varphi(n,1))^(aExp-1)*sin(varphi(n,1))*cos(theta(n,1))^bExp ...
                + 1i*ASD_theta^2*bExp*cos(theta(n,1))^(bExp-1)*sin(theta(n,1))*...
                (cos(varphi(n,1))^aExp*C-aExp*cos(varphi(n,1))^(aExp-1)*sin(varphi(n,1))*D);
            RgP(row,column) = RgP(row,column) + 1/Ncluster*A*sqrt(sigma2tilde)/ASD_varphi*...
                exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
                (X-Y*(1i*B*sigma2tilde-C*D*ASD_theta^2*sigma2tilde))*...
                exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
                exp(-B^2*sigma2tilde/2)/scaleeB(1);
        end
    end
end

RisoBS = (RisoBS+RisoBS')+eye(M);
RgP = (RgP+RgP')+eye(M);

%% Generate spatial correlation matrices for RIS
for row = 1:N
    row
    for column = row+1:N

        distanceY = NposY(row)-NposY(column);
        distanceZ = NposZ(row)-NposZ(column);
        RisoRIS(row,column) =  sin(2*pi*sqrt(distanceY^2+distanceZ^2))/(2*pi*sqrt(distanceY^2+distanceZ^2));


        for n = 1:Ncluster


            A = exp(1i*2*pi*(distanceY*sin(varphi(n,2))*cos(theta(n,2))+distanceZ*sin(theta(n,2))));
            B = 2*pi*distanceY*cos(varphi(n,2))*cos(theta(n,2));
            C = -2*pi*distanceY*cos(varphi(n,2))*sin(theta(n,2));
            D = -2*pi*distanceY*sin(varphi(n,2))*sin(theta(n,2))+2*pi*distanceZ*cos(theta(n,2));

            sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
            X = cos(varphi(n,2))^aExp*(cos(theta(n,2))^bExp-1i*bExp*cos(theta(n,2))^(bExp-1)*sin(theta(n,2))*ASD_theta^2*D);
            Y = aExp*cos(varphi(n,2))^(aExp-1)*sin(varphi(n,2))*cos(theta(n,2))^bExp ...
                + 1i*ASD_theta^2*bExp*cos(theta(n,2))^(bExp-1)*sin(theta(n,2))*...
                (cos(varphi(n,2))^aExp*C-aExp*cos(varphi(n,2))^(aExp-1)*sin(varphi(n,2))*D);
            Rh(row,column) = Rh(row,column) + 1/Ncluster*A*sqrt(sigma2tilde)/ASD_varphi*...
                exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
                (X-Y*(1i*B*sigma2tilde-C*D*ASD_theta^2*sigma2tilde))*...
                exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
                exp(-B^2*sigma2tilde/2)/scaleeB(2);


            A = exp(1i*2*pi*(distanceY*sin(varphi(n,3))*cos(theta(n,3))+distanceZ*sin(theta(n,3))));
            B = 2*pi*distanceY*cos(varphi(n,3))*cos(theta(n,3));
            C = -2*pi*distanceY*cos(varphi(n,3))*sin(theta(n,3));
            D = -2*pi*distanceY*sin(varphi(n,3))*sin(theta(n,3))+2*pi*distanceZ*cos(theta(n,3));

            sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
            X = cos(varphi(n,3))^aExp*(cos(theta(n,3))^bExp-1i*bExp*cos(theta(n,3))^(bExp-1)*sin(theta(n,3))*ASD_theta^2*D);
            Y = aExp*cos(varphi(n,3))^(aExp-1)*sin(varphi(n,3))*cos(theta(n,3))^bExp ...
                + 1i*ASD_theta^2*bExp*cos(theta(n,3))^(bExp-1)*sin(theta(n,3))*...
                (cos(varphi(n,3))^aExp*C-aExp*cos(varphi(n,3))^(aExp-1)*sin(varphi(n,3))*D);
            Rg(row,column) = Rg(row,column) + 1/Ncluster*A*sqrt(sigma2tilde)/ASD_varphi*...
                exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
                (X-Y*(1i*B*sigma2tilde-C*D*ASD_theta^2*sigma2tilde))*...
                exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
                exp(-B^2*sigma2tilde/2)/scaleeB(3);


            A = exp(1i*2*pi*(distanceY*sin(varphi(n,4))*cos(theta(n,4))+distanceZ*sin(theta(n,4))));
            B = 2*pi*distanceY*cos(varphi(n,4))*cos(theta(n,4));
            C = -2*pi*distanceY*cos(varphi(n,4))*sin(theta(n,4));
            D = -2*pi*distanceY*sin(varphi(n,4))*sin(theta(n,4))+2*pi*distanceZ*cos(theta(n,4));

            sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
            X = cos(varphi(n,4))^aExp*(cos(theta(n,4))^bExp-1i*bExp*cos(theta(n,4))^(bExp-1)*sin(theta(n,4))*ASD_theta^2*D);
            Y = aExp*cos(varphi(n,4))^(aExp-1)*sin(varphi(n,4))*cos(theta(n,4))^bExp ...
                + 1i*ASD_theta^2*bExp*cos(theta(n,4))^(bExp-1)*sin(theta(n,4))*...
                (cos(varphi(n,4))^aExp*C-aExp*cos(varphi(n,4))^(aExp-1)*sin(varphi(n,4))*D);
            Rw(row,column) = Rw(row,column) + 1/Ncluster*A*sqrt(sigma2tilde)/ASD_varphi*...
                exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
                (X-Y*(1i*B*sigma2tilde-C*D*ASD_theta^2*sigma2tilde))*...
                exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
                exp(-B^2*sigma2tilde/2)/scaleeB(4);

        end
    end
end

RisoRIS = (RisoRIS+RisoRIS')+eye(N);
Rh = (Rh+Rh')+eye(N);
Rg = (Rg+Rg')+eye(N);
Rw = (Rw+Rw')+eye(N);

%% Eigenvalue Decomposition

[UisoBS, DisoBS] = eig(RisoBS);
[UisoRIS, DisoRIS] = eig(RisoRIS.*RisoRIS);

[UgP, DgP] = eig(RgP);
[Uhg, Dhg] = eig(Rh.*Rg);
[Uwg, Dwg] = eig(Rw.*Rg);


[Uh, Dh] = eig(Rh);
[Ug, Dg] = eig(Rg);
[Uw, Dw] = eig(Rw);

%To prevent any numerical issues, equate very small negative eigenvalues
%to zero
DisoBS(DisoBS<0) = 0;
DisoRIS(DisoRIS<0) = 0;
DgP(DgP<0) = 0;
Dhg(Dhg<0) = 0;
Dwg(Dwg<0) = 0;
Dh(Dh<0) = 0;
Dg(Dg<0) = 0;
Dw(Dw<0) = 0;

RhSqrt = sqrtm(Uh*Dh*Uh');
RgSqrt = sqrtm(Ug*Dg*Ug');
RwSqrt = sqrtm(Uw*Dw*Uw');

%Determine the ranks
dd = sort(diag(DisoBS),'descend');
for rr = 1:M
    if sum(dd(1:rr))/sum(dd)>0.999999
        rank_isoBS = rr;
        break
    end
end


dd = sort(diag(DisoRIS),'descend');
for rr = 1:N
    if sum(dd(1:rr))/sum(dd)>0.999999
        rank_isoRIS = rr;
        break
    end
end

dd = sort(diag(DgP),'descend');
for rr = 1:M
    if sum(dd(1:rr))/sum(dd)>0.999999
        rank_gP = rr;
        break
    end
end


dd = sort(diag(Dhg),'descend');
for rr = 1:N
    if sum(dd(1:rr))/sum(dd)>0.999999
        rank_hg = rr;
        break
    end
end


dd = sort(diag(Dwg),'descend');
for rr = 1:N
    if sum(dd(1:rr))/sum(dd)>0.999999
        rank_wg = rr;
        break
    end
end


taup =  N/2; %pilot length

[ddgP, indexx] = sort(diag(DgP),'descend');
UgP = UgP(:,indexx);
DgP = diag(ddgP);

[ddhg, indexx] = sort(diag(Dhg),'descend');
Uhg = Uhg(:,indexx);
Dhg = diag(ddhg);

[ddwg, indexx] = sort(diag(Dwg),'descend');
Uwg = Uwg(:,indexx);
Dwg = diag(ddwg);

Uhg1 = Uhg(:,1:taup);
Uhg1b = Uhg(:,1:rank_hg);

[ddisoBS, indexx] = sort(diag(DisoBS),'descend');
UisoBS = UisoBS(:,indexx);
DisoBS = diag(ddisoBS);

[ddisoRIS, indexx] = sort(diag(DisoRIS),'descend');
UisoRIS = UisoRIS(:,indexx);
DisoRIS = diag(ddisoRIS);

UisoBS1 = UisoBS(:,1:rank_isoBS);

UisoRIS1 = UisoRIS(:,1:rank_isoRIS);

%% Channel estimation
%Carrier frequency
fc = 3e9;

%Speed of light
c = 3e8;

%Wavelength
lambda = c/fc;

distance1 = 10;
distance2 = 40;

N0 = db2pow(-174)/1000;
spacing = lambda/8;
lossComparedToIsotropic = spacing^2/(lambda^2/(4*pi));

Bandwidth = 10^5;
pathloss1 = lossComparedToIsotropic*db2pow(-35-30*log10(distance1));
pathloss2 = lossComparedToIsotropic*db2pow(-35-30*log10(distance2));

P = 0.2;
SNR = pow2db(P*pathloss1*pathloss2/(Bandwidth*N0));

SIR = -5:5:15;


nmseLMMSE_random = zeros(length(SIR),length(SNR),2);
nmseLMMSE_dft = zeros(length(SIR),length(SNR),2);
nmseLMMSE_optimized_noEMI = zeros(length(SIR),length(SNR),3);
nmseLMMSE_optimized_EMI = zeros(length(SIR),length(SNR),3);
nmseLMMSE_bound_EMI = zeros(length(SIR),length(SNR),2);
nmseLMMSE_optimized_EMI_quantized1 = zeros(length(SIR),length(SNR),2);
nmseLMMSE_optimized_EMI_quantized2 = zeros(length(SIR),length(SNR),2);


nmseLS_dft = zeros(length(SIR),length(SNR),2);
nmseLS_random = zeros(length(SIR),length(SNR),2);


nmseRSLS_iso_random = zeros(length(SIR),length(SNR),3);
nmseRSLS_iso_dft = zeros(length(SIR),length(SNR),2);
nmseRSLS_iso_optimized = zeros(length(SIR),length(SNR),3);
nmseRSLS_iso_MM = zeros(length(SIR),length(SNR),3);
nmseRSLS_iso_MM_quantized1 = zeros(length(SIR),length(SNR),2);
nmseRSLS_iso_MM_quantized2 = zeros(length(SIR),length(SNR),2);
nmseRSLS_iso_bound = zeros(length(SIR),length(SNR),2);

nmseCRLB = zeros(length(SIR),length(SNR));

Phase_dft = dftmtx(N);
Phase_dft = Phase_dft(1:taup,:);
Phase_random = exp(1i*2*pi*rand(taup,N));

DFTtaup = dftmtx(taup)/sqrt(taup);


Phase_RSLS_iso = sqrt(N*taup/rank_isoRIS)*DFTtaup(:,1:rank_isoRIS)*UisoRIS1';

Phase_RSLS_isoy = Phase_RSLS_iso;
for iterr = 1:10
    tempory = inv(UisoRIS1'*Phase_RSLS_isoy'*Phase_RSLS_isoy*UisoRIS1);
    tempory2 = Phase_RSLS_isoy*UisoRIS1;
    lambda1 = 3*trace(tempory)^2;
    A0 = -tempory2*tempory*tempory-lambda1*tempory2;
    ss =  A0*UisoRIS1';
    ss = ss(:);

    ss2 = tempory2*UisoRIS1';

    tempory4 = lambda1*ss2(:)+ss-lambda1*Phase_RSLS_isoy(:);
    Phase_RSLS_isoy = reshape(exp(1i*angle((-tempory4))),[taup,N]);
    trace(tempory)
end


    Phase_RSLS_isoy1 = reshape(exp(1i*angle((-tempory4))),[taup,N]);
    % Specify the number of bits for quantization
    b = 1;

    % Calculate the number of quantization levels
    numLevels = 2^b;

    % Calculate the quantization levels (phase shifts)
    quantizationLevels = (0:numLevels-1) * (2*pi/numLevels);

    % Initialize the quantized matrix
    quantizedA = zeros(size(Phase_RSLS_isoy1));

    % Loop through each element of the matrix
    for i = 1:size(Phase_RSLS_isoy1, 1)
        for j = 1:size(Phase_RSLS_isoy1, 2)
            % Extract the phase of the current element
            theta = wrapTo2Pi(angle(Phase_RSLS_isoy1(i, j)));

            % Find the nearest quantization level
            [~, idx] = min(abs(theta - quantizationLevels));

            % Assign the quantized phase to the new matrix
            quantizedA(i, j) = exp(1i * quantizationLevels(idx));
        end
    end

    Phase_RSLS_isoy1 = quantizedA;
    




    Phase_RSLS_isoy2 = reshape(exp(1i*angle((-tempory4))),[taup,N]);
    % Specify the number of bits for quantization
    b = 2;

    % Calculate the number of quantization levels
    numLevels = 2^b;

    % Calculate the quantization levels (phase shifts)
    quantizationLevels = (0:numLevels-1) * (2*pi/numLevels);

    % Initialize the quantized matrix
    quantizedA = zeros(size(Phase_RSLS_isoy2));

    % Loop through each element of the matrix
    for i = 1:size(Phase_RSLS_isoy2, 1)
        for j = 1:size(Phase_RSLS_isoy2, 2)
            % Extract the phase of the current element
            theta = wrapTo2Pi(angle(Phase_RSLS_isoy2(i, j)));

            % Find the nearest quantization level
            [~, idx] = min(abs(theta - quantizationLevels));

            % Assign the quantized phase to the new matrix
            quantizedA(i, j) = exp(1i * quantizationLevels(idx));
        end
    end

    Phase_RSLS_isoy2 = quantizedA;


Phase_RSLS_isox = exp(1i*angle(Phase_RSLS_iso));
tempory = inv(UisoRIS1'*Phase_RSLS_isox'*Phase_RSLS_isox*UisoRIS1);
trace(tempory)
RgP = UgP*DgP*UgP';
Rhg = Uhg*Dhg*Uhg';
Rwg = Uwg*Dwg*Uwg';

RgPSqrt = sqrtm(RgP);

R_BSiso = UisoBS*DisoBS*UisoBS';
R_RISiso = UisoRIS*DisoRIS*UisoRIS';

BLMMSE_random = kron(RgP, Phase_random*Rhg);
BLMMSE_dft = kron(RgP, Phase_dft*Rhg);

temporb = (Phase_random'*Phase_random);
aLS_random = real(M*trace(eye(size(temporb,1))/temporb));
mLS_random = kron(eye(M),temporb\Phase_random');



temporb = (Phase_dft'*Phase_dft);
aLS_dft = real(M*trace(eye(size(temporb,1))/temporb));
mLS_dft = kron(eye(M),temporb\Phase_dft');

bLS_random = trace(RgP)*trace(Rwg);
bLS_dft = bLS_random;


tempora = UisoRIS1'*(Phase_random'*Phase_random);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_random = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_random = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_random = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_random')));


tempora = UisoRIS1'*(Phase_dft'*Phase_dft);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_dft = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_dft = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_dft = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_dft')));


tempora = UisoRIS1'*(Phase_RSLS_isox'*Phase_RSLS_isox);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_optimized = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_optimized = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_optimized = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_RSLS_isox')));


tempora = UisoRIS1'*(Phase_RSLS_isoy'*Phase_RSLS_isoy);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_MM = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_MM = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_MM = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_RSLS_isoy')));

tempora = UisoRIS1'*(Phase_RSLS_isoy1'*Phase_RSLS_isoy1);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_MM_quantized1 = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_MM_quantized1 = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_MM_quantized1 = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_RSLS_isoy1')));

tempora = UisoRIS1'*(Phase_RSLS_isoy2'*Phase_RSLS_isoy2);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_MM_quantized2 = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_MM_quantized2 = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_MM_quantized2 = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_RSLS_isoy2')));

tempora = UisoRIS1'*(Phase_RSLS_iso'*Phase_RSLS_iso);
temporb = (tempora*UisoRIS1);
temporc = eye(size(temporb,1))/temporb;
aRSLS_iso_bound = real(rank_isoBS*trace(temporc));
tempord = temporb\tempora;
bRSLS_iso_bound = real(rank_isoBS*trace(tempord*Rwg*tempord'));
mRSLS_iso_bound = kron(UisoBS1*UisoBS1',UisoRIS1*(temporb\(UisoRIS1'*Phase_RSLS_iso')));



for scenSIR = 1:length(SIR)
    sir0 = 10^(SIR(scenSIR)/10);


    Bdot = (Rhg + 1/sir0*Rwg);

    [UBdot, DBdot] = eig(Bdot);
    DBdot = real(DBdot);
    DBdot(DBdot<0) = 0;
    dd = sort(diag(DBdot),'descend');
    for rr = 1:N
        if sum(dd(1:rr))/sum(dd)>0.999999
            rank_Bdot = rr;
            break
        end
    end

    [ddBdot, indexx] = sort(diag(DBdot),'descend');
    UBdot = UBdot(:,indexx);
    DBdot = diag(ddBdot);
    UBdot1 = UBdot(:,1:rank_Bdot);
    DBdot1 = diag(ddBdot(1:rank_Bdot));
    DBdot1inv05 = diag((ddBdot(1:rank_Bdot)).^(-0.5));
    Gdot = DBdot1inv05*UBdot1'*Rhg*Rhg*UBdot1*DBdot1inv05;

    Bdot = UBdot*DBdot*UBdot';
    ALMMSE_random = kron(RgP, Phase_random*Bdot*Phase_random');
    ALMMSE_dft = kron(RgP, Phase_dft*Bdot*Phase_dft');

    for scenSNR = 1:length(SNR)
        snr0 = 10^(SNR(scenSNR)/10);
        [sir0 snr0]
        cvx_begin quiet
        variable lambda2(taup,1)
        variable sss(rank_Bdot,1)
        minimize sum(sss)
        subject to

        for tt = 1:rank_Bdot
            sum(ddgP*real(Gdot(tt,tt)).*inv_pos(snr0*ddgP*ddBdot(tt)*lambda2(tt)+1))<=sss(tt);
        end

        lambda2>=zeros(taup,1);
        sum(lambda2)<=N*taup;
        cvx_end

        Phase_LMMSE_EMI = dftmtx(taup)/sqrt(taup)*sqrt(diag(lambda2))*UBdot(:,1:taup)';

        Phase_LMMSEx_EMI = exp(1i*angle(Phase_LMMSE_EMI));

        % Specify the number of bits for quantization
        b = 1;

        % Calculate the number of quantization levels
        numLevels = 2^b;

        % Calculate the quantization levels (phase shifts)
        quantizationLevels = (0:numLevels-1) * (2*pi/numLevels);

        % Initialize the quantized matrix
        quantizedA = zeros(size(Phase_LMMSEx_EMI));

        % Loop through each element of the matrix
        for i = 1:size(Phase_LMMSEx_EMI, 1)
            for j = 1:size(Phase_LMMSEx_EMI, 2)
                % Extract the phase of the current element
                theta = angle(Phase_LMMSEx_EMI(i, j));

                % Find the nearest quantization level
                [~, idx] = min(abs(theta - quantizationLevels));

                % Assign the quantized phase to the new matrix
                quantizedA(i, j) = exp(1i * quantizationLevels(idx));
            end
        end

        Phase_LMMSEx_EMI1 = quantizedA;


        % Specify the number of bits for quantization
        b = 2;

        % Calculate the number of quantization levels
        numLevels = 2^b;

        % Calculate the quantization levels (phase shifts)
        quantizationLevels = (0:numLevels-1) * (2*pi/numLevels);

        % Initialize the quantized matrix
        quantizedA = zeros(size(Phase_LMMSEx_EMI));

        % Loop through each element of the matrix
        for i = 1:size(Phase_LMMSEx_EMI, 1)
            for j = 1:size(Phase_LMMSEx_EMI, 2)
                % Extract the phase of the current element
                theta = angle(Phase_LMMSEx_EMI(i, j));

                % Find the nearest quantization level
                [~, idx] = min(abs(theta - quantizationLevels));

                % Assign the quantized phase to the new matrix
                quantizedA(i, j) = exp(1i * quantizationLevels(idx));
            end
        end

        Phase_LMMSEx_EMI2 = quantizedA;

        
        %%%%%%%%%%%%%%%%%

        cvx_begin quiet
        variable lambda2(taup,1)
        variable sss(taup,1)
        minimize sum(sss)
        subject to

        for tt = 1:taup
            sum(ddgP*ddhg(tt).*inv_pos(snr0*ddgP*ddhg(tt)*lambda2(tt)+1))<=sss(tt);
        end

        lambda2>=zeros(taup,1);
        sum(lambda2)<=N*taup;
        cvx_end

        Phase_LMMSE_noEMI = dftmtx(taup)/sqrt(taup)*sqrt(diag(lambda2))*Uhg1';

        Phase_LMMSEx_noEMI = exp(1i*angle(Phase_LMMSE_noEMI));

        %%%%%%%%%%%%%%%%%%%%%555
        ALMMSE_optimized_noEMI = kron(RgP, Phase_LMMSEx_noEMI*Bdot*Phase_LMMSEx_noEMI');
        BLMMSE_optimized_noEMI = kron(RgP, Phase_LMMSEx_noEMI*Rhg);
        temporCRB2= kron(eye(M),Phase_RSLS_isoy);
        temporCRB = temporCRB2*kron(UisoBS1,UisoRIS1);
        CRLB = trace(inv(snr0*temporCRB'*((temporCRB2*snr0/sir0*kron(RgP,Rwg)*temporCRB2'+eye(M*taup))\temporCRB)));

        ALMMSE_optimized_EMI = kron(RgP, Phase_LMMSEx_EMI*Bdot*Phase_LMMSEx_EMI');
        ALMMSE_optimized_EMI_quantized1 = kron(RgP, Phase_LMMSEx_EMI1*Bdot*Phase_LMMSEx_EMI1');
        ALMMSE_optimized_EMI_quantized2 = kron(RgP, Phase_LMMSEx_EMI2*Bdot*Phase_LMMSEx_EMI2');

        ALMMSE_bound_EMI = kron(RgP, Phase_LMMSE_EMI*Bdot*Phase_LMMSE_EMI');

        BLMMSE_optimized_EMI = kron(RgP, Phase_LMMSEx_EMI*Rhg);
        BLMMSE_optimized_EMI_quantized1 = kron(RgP, Phase_LMMSEx_EMI1*Rhg);
        BLMMSE_optimized_EMI_quantized2 = kron(RgP, Phase_LMMSEx_EMI2*Rhg);

        BLMMSE_bound_EMI = kron(RgP, Phase_LMMSE_EMI*Rhg);


        nmseLMMSE_random(scenSIR,scenSNR,1) = real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_random'*((ALMMSE_random*snr0+eye(M*taup))\BLMMSE_random))/M/N);
        nmseLMMSE_dft(scenSIR,scenSNR,1) = real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_dft'*((ALMMSE_dft*snr0+eye(M*taup))\BLMMSE_dft))/M/N);
        nmseLMMSE_optimized_noEMI(scenSIR,scenSNR,1) = real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_optimized_noEMI'*((ALMMSE_optimized_noEMI*snr0+eye(M*taup))\BLMMSE_optimized_noEMI))/M/N);

        nmseLMMSE_optimized_EMI(scenSIR,scenSNR,1) = real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_optimized_EMI'*((ALMMSE_optimized_EMI*snr0+eye(M*taup))\BLMMSE_optimized_EMI))/M/N);
        nmseLMMSE_optimized_EMI_quantized1(scenSIR,scenSNR,1) = real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_optimized_EMI_quantized1'*((ALMMSE_optimized_EMI_quantized1*snr0+eye(M*taup))\BLMMSE_optimized_EMI_quantized1))/M/N);
        nmseLMMSE_optimized_EMI_quantized2(scenSIR,scenSNR,1) = real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_optimized_EMI_quantized2'*((ALMMSE_optimized_EMI_quantized2*snr0+eye(M*taup))\BLMMSE_optimized_EMI_quantized2))/M/N);

        nmseLMMSE_bound_EMI(scenSIR,scenSNR,1) =  real(trace(RgP)*trace(Rhg)/M/N-snr0*trace(BLMMSE_bound_EMI'*((ALMMSE_bound_EMI*snr0+eye(M*taup))\BLMMSE_bound_EMI))/M/N);

        nmseLS_dft(scenSIR,scenSNR,1) = aLS_dft/snr0/M/N + bLS_dft/sir0/M/N;


        nmseRSLS_iso_dft(scenSIR,scenSNR,1) = aRSLS_iso_dft/snr0/M/N + bRSLS_iso_dft/sir0/M/N ;
        nmseRSLS_iso_optimized(scenSIR,scenSNR,1) = aRSLS_iso_optimized/snr0/M/N + bRSLS_iso_optimized/sir0/M/N;
        nmseRSLS_iso_MM(scenSIR,scenSNR,1) = aRSLS_iso_MM/snr0/M/N + bRSLS_iso_MM/sir0/M/N;
        nmseRSLS_iso_MM_quantized1(scenSIR,scenSNR,1) = aRSLS_iso_MM_quantized1/snr0/M/N + bRSLS_iso_MM_quantized1/sir0/M/N;
        nmseRSLS_iso_MM_quantized2(scenSIR,scenSNR,1) = aRSLS_iso_MM_quantized2/snr0/M/N + bRSLS_iso_MM_quantized2/sir0/M/N;

        nmseRSLS_iso_bound(scenSIR,scenSNR,1) = aRSLS_iso_bound/snr0/M/N + bRSLS_iso_bound/sir0/M/N;
        nmseCRLB(scenSIR,scenSNR) =  CRLB/(M*N);
        nbrOfRealizations = 100;

        for channelTrial = 1:nbrOfRealizations
            channelTrial
            channel_h = RhSqrt*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            EMI_w = RwSqrt*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));
            EMI_w2 = RwSqrt*sqrt(0.5)*(randn(N,taup)+1i*randn(N,taup));
            channel_G = RgSqrt*sqrt(0.5)*(randn(N,M)+1i*randn(N,M))*conj(RgPSqrt);
            channel_x = zeros(M*N,1);
            EMI_x = zeros(M*N,1);
            EMI_x2 = zeros(M*N,taup);
            receivedSignal_LMMSE_random = zeros(M*taup,1);
            receivedSignal_LMMSE_dft = zeros(M*taup,1);
            receivedSignal_LMMSE_optimized_noEMI = zeros(M*taup,1);
            receivedSignal_LMMSE_optimized_noEMI2 = zeros(M*taup,1);

            receivedSignal_LMMSE_optimized_EMI = zeros(M*taup,1);
            receivedSignal_LMMSE_optimized_EMI_quantized1  = zeros(M*taup,1);
            receivedSignal_LMMSE_optimized_EMI_quantized2  = zeros(M*taup,1);

            receivedSignal_LMMSE_optimized_EMI2 = zeros(M*taup,1);

            receivedSignal_LMMSE_bound_EMI = zeros(M*taup,1);
            receivedSignal_LS_random = zeros(M*taup,1);
            receivedSignal_LS_dft = zeros(M*taup,1);
            receivedSignal_RSLS_iso_random = zeros(M*taup,1);
            receivedSignal_RSLS_iso_random2 = zeros(M*taup,1);
            receivedSignal_RSLS_iso_dft = zeros(M*taup,1);

            receivedSignal_RSLS_iso_optimized = zeros(M*taup,1);
            receivedSignal_RSLS_iso_optimized2 = zeros(M*taup,1);

            receivedSignal_RSLS_iso_MM = zeros(M*taup,1);
            receivedSignal_RSLS_iso_MM2 = zeros(M*taup,1);
            receivedSignal_RSLS_iso_MM_quantized1 = zeros(M*taup,1);
            receivedSignal_RSLS_iso_MM_quantized2 = zeros(M*taup,1);

            receivedSignal_RSLS_iso_bound = zeros(M*taup,1);

            for mm = 1:M
                channel_x((mm-1)*N+1:mm*N,1) = channel_h.*channel_G(:,mm);
                EMI_x((mm-1)*N+1:mm*N,1) = EMI_w.*channel_G(:,mm);
                EMI_x2((mm-1)*N+1:mm*N,:) = EMI_w2.*channel_G(:,mm);

                receivedSignal_LMMSE_random((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_random*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_random*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LMMSE_dft((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_dft*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_dft*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LMMSE_optimized_noEMI((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSEx_noEMI*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_LMMSEx_noEMI*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LMMSE_optimized_noEMI2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSEx_noEMI*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*sum(Phase_LMMSEx_noEMI.*EMI_x2((mm-1)*N+1:mm*N,:).',2) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));



                receivedSignal_LMMSE_optimized_EMI((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSEx_EMI*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_LMMSEx_EMI*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LMMSE_optimized_EMI_quantized1((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSEx_EMI1*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_LMMSEx_EMI1*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));
                receivedSignal_LMMSE_optimized_EMI_quantized2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSEx_EMI2*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_LMMSEx_EMI2*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1)); 
               

                receivedSignal_LMMSE_optimized_EMI2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSEx_EMI*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*sum(Phase_LMMSEx_EMI.*EMI_x2((mm-1)*N+1:mm*N,:).',2) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LMMSE_bound_EMI((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_LMMSE_EMI*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_LMMSE_EMI*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LS_random((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_random*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_random*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_LS_dft((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_dft*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_dft*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_RSLS_iso_random((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_random*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_random*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_RSLS_iso_random2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_random*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*sum(Phase_random.*EMI_x2((mm-1)*N+1:mm*N,1).',2) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));


                receivedSignal_RSLS_iso_dft((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_dft*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_dft*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));


                receivedSignal_RSLS_iso_optimized((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_isox*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_RSLS_isox*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_RSLS_iso_optimized2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_isox*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*sum(Phase_RSLS_isox.*EMI_x2((mm-1)*N+1:mm*N,:).',2) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));


                receivedSignal_RSLS_iso_MM((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_isoy*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_RSLS_isoy*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_RSLS_iso_MM2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_isoy*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*sum(Phase_RSLS_isoy.*EMI_x2((mm-1)*N+1:mm*N,:).',2) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

               receivedSignal_RSLS_iso_MM_quantized1((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_isoy1*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_RSLS_isoy1*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                 receivedSignal_RSLS_iso_MM_quantized2((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_isoy2*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_RSLS_isoy2*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

                receivedSignal_RSLS_iso_bound((mm-1)*taup+1:mm*taup,1) = ...
                    sqrt(snr0)*Phase_RSLS_iso*channel_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(snr0/sir0)*Phase_RSLS_iso*EMI_x((mm-1)*N+1:mm*N,1) ...
                    +sqrt(0.5)*(randn(taup,1)+1i*randn(taup,1));

            end


            estimate_LMMSE_random = sqrt(snr0)*BLMMSE_random'*((ALMMSE_random*snr0+eye(M*taup))\receivedSignal_LMMSE_random);
            estimate_LMMSE_dft = sqrt(snr0)*BLMMSE_dft'*((ALMMSE_dft*snr0+eye(M*taup))\receivedSignal_LMMSE_dft);
            estimate_LMMSE_optimized_noEMI = sqrt(snr0)*BLMMSE_optimized_noEMI'*((ALMMSE_optimized_noEMI*snr0+eye(M*taup))\receivedSignal_LMMSE_optimized_noEMI);
            estimate_LMMSE_optimized_noEMI2 = sqrt(snr0)*BLMMSE_optimized_noEMI'*((ALMMSE_optimized_noEMI*snr0+eye(M*taup))\receivedSignal_LMMSE_optimized_noEMI2);

            estimate_LMMSE_optimized_EMI = sqrt(snr0)*BLMMSE_optimized_EMI'*((ALMMSE_optimized_EMI*snr0+eye(M*taup))\receivedSignal_LMMSE_optimized_EMI);
            estimate_LMMSE_optimized_EMI2 = sqrt(snr0)*BLMMSE_optimized_EMI'*((ALMMSE_optimized_EMI*snr0+eye(M*taup))\receivedSignal_LMMSE_optimized_EMI2);
            estimate_LMMSE_optimized_EMI_quantized1 = sqrt(snr0)*BLMMSE_optimized_EMI_quantized1'*((ALMMSE_optimized_EMI_quantized1*snr0+eye(M*taup))\receivedSignal_LMMSE_optimized_EMI_quantized1);
            estimate_LMMSE_optimized_EMI_quantized2 = sqrt(snr0)*BLMMSE_optimized_EMI_quantized2'*((ALMMSE_optimized_EMI_quantized2*snr0+eye(M*taup))\receivedSignal_LMMSE_optimized_EMI_quantized2);

            estimate_LMMSE_bound_EMI = sqrt(snr0)*BLMMSE_bound_EMI'*((ALMMSE_bound_EMI*snr0+eye(M*taup))\receivedSignal_LMMSE_bound_EMI);
            estimate_LS_random = mLS_random*receivedSignal_LS_random/sqrt(snr0);
            estimate_LS_dft = mLS_dft*receivedSignal_LS_dft/sqrt(snr0);
            estimate_RSLS_iso_random = mRSLS_iso_random*receivedSignal_RSLS_iso_random/sqrt(snr0);
            estimate_RSLS_iso_random2 = mRSLS_iso_random*receivedSignal_RSLS_iso_random2/sqrt(snr0);

            estimate_RSLS_iso_dft = mRSLS_iso_dft*receivedSignal_RSLS_iso_dft/sqrt(snr0);

            estimate_RSLS_iso_optimized = mRSLS_iso_optimized*receivedSignal_RSLS_iso_optimized/sqrt(snr0);
            estimate_RSLS_iso_optimized2 = mRSLS_iso_optimized*receivedSignal_RSLS_iso_optimized2/sqrt(snr0);
            estimate_RSLS_iso_MM = mRSLS_iso_MM*receivedSignal_RSLS_iso_MM/sqrt(snr0);
            estimate_RSLS_iso_MM2 = mRSLS_iso_MM*receivedSignal_RSLS_iso_MM2/sqrt(snr0);
            estimate_RSLS_iso_bound = mRSLS_iso_bound*receivedSignal_RSLS_iso_bound/sqrt(snr0);
            estimate_RSLS_iso_MM_quantized1 = mRSLS_iso_MM_quantized1*receivedSignal_RSLS_iso_MM_quantized1/sqrt(snr0);
            estimate_RSLS_iso_MM_quantized2 = mRSLS_iso_MM_quantized2*receivedSignal_RSLS_iso_MM_quantized2/sqrt(snr0);


            nmseLMMSE_random(scenSIR,scenSNR,2) = nmseLMMSE_random(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_random)^2/M/N/nbrOfRealizations;
            nmseLMMSE_dft(scenSIR,scenSNR,2) = nmseLMMSE_dft(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_dft)^2/M/N/nbrOfRealizations;
            nmseLMMSE_optimized_noEMI(scenSIR,scenSNR,2) = nmseLMMSE_optimized_noEMI(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_optimized_noEMI)^2/M/N/nbrOfRealizations;
            nmseLMMSE_optimized_noEMI(scenSIR,scenSNR,3) = nmseLMMSE_optimized_noEMI(scenSIR,scenSNR,3) + norm(channel_x-estimate_LMMSE_optimized_noEMI2)^2/M/N/nbrOfRealizations;

            nmseLMMSE_optimized_EMI(scenSIR,scenSNR,2) = nmseLMMSE_optimized_EMI(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_optimized_EMI)^2/M/N/nbrOfRealizations;
            nmseLMMSE_optimized_EMI_quantized1(scenSIR,scenSNR,2) = nmseLMMSE_optimized_EMI_quantized1(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_optimized_EMI_quantized1)^2/M/N/nbrOfRealizations;
            nmseLMMSE_optimized_EMI_quantized2(scenSIR,scenSNR,2) = nmseLMMSE_optimized_EMI_quantized2(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_optimized_EMI_quantized2)^2/M/N/nbrOfRealizations;


            nmseLMMSE_optimized_EMI(scenSIR,scenSNR,3) = nmseLMMSE_optimized_EMI(scenSIR,scenSNR,3) + norm(channel_x-estimate_LMMSE_optimized_EMI2)^2/M/N/nbrOfRealizations;

            nmseLMMSE_bound_EMI(scenSIR,scenSNR,2) = nmseLMMSE_bound_EMI(scenSIR,scenSNR,2) + norm(channel_x-estimate_LMMSE_bound_EMI)^2/M/N/nbrOfRealizations;
            nmseLS_random(scenSIR,scenSNR,2) = nmseLS_random(scenSIR,scenSNR,2) + norm(channel_x-estimate_LS_random)^2/M/N/nbrOfRealizations;

            nmseLS_dft(scenSIR,scenSNR,2) = nmseLS_dft(scenSIR,scenSNR,2) + norm(channel_x-estimate_LS_dft)^2/M/N/nbrOfRealizations;

            nmseRSLS_iso_random(scenSIR,scenSNR,2) = nmseRSLS_iso_random(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_random)^2/M/N/nbrOfRealizations;
            nmseRSLS_iso_random(scenSIR,scenSNR,3) = nmseRSLS_iso_random(scenSIR,scenSNR,3) + norm(channel_x-estimate_RSLS_iso_random2)^2/M/N/nbrOfRealizations;


            nmseRSLS_iso_dft(scenSIR,scenSNR,2) = nmseRSLS_iso_dft(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_dft)^2/M/N/nbrOfRealizations;

            nmseRSLS_iso_optimized(scenSIR,scenSNR,2) = nmseRSLS_iso_optimized(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_optimized)^2/M/N/nbrOfRealizations;
            nmseRSLS_iso_optimized(scenSIR,scenSNR,3) = nmseRSLS_iso_optimized(scenSIR,scenSNR,3) + norm(channel_x-estimate_RSLS_iso_optimized2)^2/M/N/nbrOfRealizations;
            nmseRSLS_iso_MM(scenSIR,scenSNR,2) = nmseRSLS_iso_MM(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_MM)^2/M/N/nbrOfRealizations;
            nmseRSLS_iso_MM(scenSIR,scenSNR,3) = nmseRSLS_iso_MM(scenSIR,scenSNR,3) + norm(channel_x-estimate_RSLS_iso_MM2)^2/M/N/nbrOfRealizations;
            nmseRSLS_iso_MM_quantized1(scenSIR,scenSNR,2) = nmseRSLS_iso_MM_quantized1(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_MM_quantized1)^2/M/N/nbrOfRealizations;
            nmseRSLS_iso_MM_quantized2(scenSIR,scenSNR,2) = nmseRSLS_iso_MM_quantized2(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_MM_quantized2)^2/M/N/nbrOfRealizations;

            nmseRSLS_iso_bound(scenSIR,scenSNR,2) = nmseRSLS_iso_bound(scenSIR,scenSNR,2) + norm(channel_x-estimate_RSLS_iso_bound)^2/M/N/nbrOfRealizations;

        end

    end
end



