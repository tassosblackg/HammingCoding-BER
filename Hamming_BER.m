%% Hamming SDD & HDD compare ber for m=3 & m=4
  %|Hamming test BER for HDD and SDD-- Implementation of Hamming (7,4) and (15,11)
  %|Implementation made by using Octave 4.2.2 on Ubuntu 16.04

close all;clear;clc;

%load pkg Ocatve -- comment out in case of Matlab
pkg load communications;
pkg load signal;
%------------------------------------------------

%%|Define Var stage|------------------------
N=10000; %number of interations monte carlo
n1=7;k1=4;  %Hammin (7,4)
n2=15;k2=11; %Hamming (15,11)

%----------------------|Hamming(7,4)|---------------------------------------------------------------------------------
[h,g]=hammgen(n1-k1); %create hamming parity and generator matrix
dictII=mod(de2bi(0:(2^k1)-1)*g,2); %create dictionary with all possible messages and their parity bits matrix(2^k1,n1)
%----------------------|Hamming(15,11)|-------------------------------------------------------------------------------
[h2,g2]=hammgen(n2-k2);
dictIII=mod(de2bi(0:(2^k2)-1)*g2,2);%create dictionary with all possible messages and their parity bits matrix(2^k2,n2)

% Create random bit_streams messages matrix N streams-messages of size k
msg4=(sign(randn(N,k1))+1)/2;
msg11=(sign(randn(N,k2))+1)/2;
msg11_m=msg11; %moduated messages
msg11_m(msg11_m==0)=-1;
% create hamming coded messages
coded_w7=mod(msg4*g,2);
coded_w7(coded_w7==0)=-1; % modulation one bit per symbol
coded_w15=mod(msg11*g2,2);
coded_w15(coded_w15==0)=-1; %modulation
%--|SNR vectors
SNRdB=[0:1:12]; %SNR vector represnted in dB
SNR=10.^(SNRdB./10);

%---|BER vectors|
ber_noCoding=zeros(1,length(SNR));
ber_7sdd=zeros(1,length(SNR));
ber_7hdd=zeros(1,length(SNR));
ber_15sdd=zeros(1,length(SNR));
ber_15hdd=zeros(1,length(SNR));

%% MONTE CARLO EXPIRMENTS
for i=1:length(SNR)
  yI=(sqrt(SNR(i))*msg11_m)+randn(N,k2);
  yII=(sqrt(SNR(i))*coded_w7)+randn(N,n1); % chanel out signal with AWGN (Nxn1) matrix
  yIII=(sqrt(SNR(i))*coded_w15)+randn(N,n2);
  for j=1:N
    %%-------- |No coding case|----------
    yI_hd=yI(j,:);
    yI_hd(yI_hd>0)=1;
    yI_hd(yI_hd<=0)=0;
    ber_noCoding(i)=ber_noCoding(i)+length(find(yI_hd~=msg11(j,:)))/k2;

    %%----------------|Hamming (7,4)|-------------------------------------------
    %SDD(7,4)
    for m=1:length(dictII)
      distanceII(m)=norm(yII(j,:)-dictII(m,:));
    end
    [minv,mini]=min(distanceII);
    decoded_w7=dictII(mini,:); %get the coded_word with minimum distance
    cw7=coded_w7(j,:);
    cw7(cw7==-1)=0; %revert modulation to calculate error between decoded and original coded message
    ber_7sdd(i)=ber_7sdd(i)+length(find(decoded_w7~=cw7))/n1;

    % HDD(7,4)
    yII_hd=yII(j,:);
    yII_hd(yII_hd>0)=1;
    yII_hd(yII_hd<=0)=0;
    zII=h*yII_hd';% calculate syndrome
    for g1=1:length(yII_hd)
      if (zII==h(:,g1))
        yII_hd(g1)=~yII_hd(g1); %fix mistake
      end
    end
    ber_7hdd(i)=ber_7hdd(i)+length(find(yII_hd~=cw7))/n1;

    %%------------------|Hamming(15,11)|----------------------------------------
    %SDD(7,4)
    for l=1:length(dictIII)
      distanceIII(l)=norm(yIII(j,:)-dictIII(l,:));
    end
    [minv2,mini2]=min(distanceIII);
    decoded_w15=dictIII(mini2,:); %get the coded_word with minimum distance
    cw15=coded_w15(j,:);
    cw15(cw15==-1)=0; %revert modulation to calculate error between decoded and original coded message
    ber_15sdd(i)=ber_15sdd(i)+length(find(decoded_w15~=cw15))/n2;

    % HDD(7,4)
    yIII_hd=yIII(j,:);
    yIII_hd(yIII_hd>0)=1;
    yIII_hd(yIII_hd<=0)=0;
    zIII=h2*yIII_hd';% calculate syndrome
    for g2=1:length(yIII_hd)
      if (zIII==h2(:,g2))
        yIII_hd(g1)=~yIII_hd(g1);
      end
    end
    ber_15hdd(i)=ber_15hdd(i)+length(find(yIII_hd~=cw15))/n2;

  end


end
%%__|Final average BER|_______
ber_noCoding=ber_noCoding/N;
ber_7sdd=ber_7sdd/N;
ber_7hdd=ber_7hdd/N;
ber_15sdd=ber_15sdd/N;
ber_15hdd=ber_15hdd/N;
%%------|PLOTS|-------------------------------------
figure(1)
semilogy(SNRdB,ber_noCoding,"b--")
hold on;
semilogy(SNRdB,ber_7sdd,"r-s")
hold on;
semilogy(SNRdB,ber_7hdd,"k-s")
hold on;
semilogy(SNRdB,ber_15sdd,"g-*")
hold on;
semilogy(SNRdB,ber_15hdd,"m-*")
title({"BER for no coding and Hamming","Hamming(7,4),Hamming(15,11)"});
ylabel({"BER"});
xlabel({"SNR(dB)"});

legend("No coding","Hamming7-SDD","Hamming7-HDD","Hamming15-SDD","Hamming15-HDD")
