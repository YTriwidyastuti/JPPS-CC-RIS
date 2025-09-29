close all
clear
clc

fprintf(datestr(datetime(now,'ConvertFrom','datenum')))
fprintf('\n')

PowSrc = db2pow(20-30);
% BW = 0.05 * (10^9);
BW = 50 * (10^6);
sigma2 = db2pow(-80-30);

M = 4;
N = 4;
U = 3;
B = 3;
kap = 1;
numdata = 1000;

load HSR_S3_R8_M4_N4.csv
inpHSR = HSR_S3_R8_M4_N4(2:N+1,:);
HSR = inpHSR(:,1:M).*exp(1i*inpHSR(:,M+1:M*2));

load HRD_RPG_S3_R8_M4_N4_U3.csv
inpHRD = HRD_RPG_S3_R8_M4_N4_U3(2:numdata*U+1,:);
allHRD = inpHRD(:,1:N).*exp(1i*inpHRD(:,N+1:N*2));

PAcc = PowSrc/U;
PSpool = wrapTo2Pi((1:2^B)*2*pi/(2^B));
numrep = 10;

SRcc = zeros([numrep,numdata]);
timeCC = zeros([numrep,numdata]);
timeSR = zeros([numrep,numdata]);

for idx = 1:numdata
hRD = transpose(allHRD(U*(idx-1)+1:U*idx,:));

for rep = 1:numrep
rng(rep)

PScc = zeros([N,1]);
tic
for nn = 1:N
    MaxSum = 0;
    for qq = 1:2^B
        Mat = real(transpose(hRD(nn,:))*kap*exp(1i*PSpool(qq))*HSR(nn,:));
        MatSum = sum(sum(Mat));
        if MaxSum <= MatSum
            MaxSum = MatSum;
            PScc(nn) = PSpool(qq);
        end
    end
end

BF = zeros([M,U]);
BF_flat = zeros([M*U,1]);
for ue = 1:U
    hRDrisHSR = transpose(hRD(:,ue))*kap*diag(exp(1i*PScc))*HSR;
    BF(:,ue) = transpose(conj(hRDrisHSR))/norm(hRDrisHSR);
    BF_flat(M*(ue-1)+1:M*ue) = BF(:,ue);
end

PAcc = zeros([U,1]);
for ue = 1:U
    PAcc(ue) = abs(transpose(hRD(:,ue))*kap*diag(exp(1i*PScc))*HSR*BF(:,ue))^2;
end
PAcc = PAcc * PowSrc / sum(PAcc);
timeCC(rep,idx) = toc;

tic
for ue = 1:U
    hRDrisHSR = transpose(hRD(:,ue))*kap*diag(exp(1i*PScc))*HSR;
    UseSig = PAcc(ue)*(abs(hRDrisHSR*BF(:,ue))^2);
    ItfSig = 0;
    for itf = 1:U
        if itf ~= ue
            ItfSig = ItfSig + PAcc(ue)*(abs(hRDrisHSR*BF(:,itf))^2);
        end
    end
    SNR = UseSig / (ItfSig + sigma2);
    SRcc(rep,idx) = SRcc(rep,idx) + BW*log2(1+SNR);
end
timeSR(rep,idx) = toc;

if mod(rep,100) == 0
    fprintf('Repeat %d ',rep)
    fprintf(datestr(datetime(now,'ConvertFrom','datenum')))
    fprintf('\n')
end

if mod(idx,100) == 0
    fprintf('Data %d ',idx)
    fprintf(datestr(datetime(now,'ConvertFrom','datenum')))
    fprintf('\n')
end

end
end

ave_time = mean(mean(timeCC));
alltime = mean(timeCC,2);
aveSRcc = mean(mean(SRcc));
allSRcc = mean(SRcc,2);

fprintf(datestr(datetime(now,'ConvertFrom','datenum')))
fprintf('\n')
