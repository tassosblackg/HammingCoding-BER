%%|Information Theory Project
  %|Test Bench to compare behaviour of hamming coding and no coding information
  %|Hamming test BER for HDD and SDD-- Implementation of Hamming (7,4) and (15,11)
  %|Implementation made by using Octave 4.2.2 on Ubuntu 16.04

close all;clear;clc;
% need this Octave---comment out this part if run in Matlab-----
pkg load communications;
pkg load signal;
%---------------------------------------------------------------
%%  |Init Part--Global Var|-----------------------------------------------------------------
bits=[1 0];
n1=7;k1=4; %Hamming (7,4) params
n2=15;k2=11; %Hamming (15,11) params
N=10000; %number of montecarlo recursions
SNRdB=[0:3:12]; %SNR vector represnted in dB
SNR=10.^(SNRdB./10);
ber_noCoding=zeros(1,length(SNR)); %only for hdd
ber_7hdd=zeros(1,length(SNR));
ber_7sdd=zeros(1,length(SNR));
ber_15hdd=zeros(1,length(SNR));
ber_15sdd=zeros(1,length(SNR));

% |Hamming (7,4)|----------------------------------------------------------------
G=[1 0 0 0 0 1 1; 0 1 0 0 1 0 1; 0 0 1 0 1 1 0; 0 0 0 1 1 1 1]; %generator matrix ,coded_word=word*G--parity bits at the end
H=[0 1 1 1 1 0 0; 1 0 1 1 0 1 0; 1 1 0 1 0 0 1];  %Parity matrix
dictII=de2bi(0:2^k1-1); %all possible messages with length k1
%add parity bits at the end
dictII(1:length(dictII),5)=xor(dictII(:,1),xor(dictII(:,2),dictII(:,4)));
dictII(1:length(dictII),6)=xor(dictII(:,1),xor(dictII(:,3),dictII(:,4)));
dictII(1:length(dictII),7)=xor(dictII(:,2),xor(dictII(:,3),dictII(:,4)));
% |Hamming (15,11)|-------------------------------------------------------------
G2=[1 0 0 0 0 0 0 0 0 0 0 1 1 0 0; 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0;
    0 0 0 1 0 0 0 0 0 0 0 1 1 1 0; 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1
    0 0 0 0 0 0 1 0 0 0 0 1 1 0 1; 0 0 0 0 0 0 0 1 0 0 0 0 0 1 1; 0 0 0 0 0 0 0 0 1 0 0 1 0 1 1
    0 0 0 0 0 0 0 0 0 1 0 0 1 1 1; 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1];
H2=[1 0 1 1 1 0 0 0 1 1 1 1 0 0 0; 1 1 0 1 1 0 1 1 0 0 1 0 1 0 0;
    1 1 1 0 1 1 0 1 1 0 0 0 0 1 0; 1 1 1 1 0 1 1 0 0 1 0 0 0 0 1];
dictIII=de2bi(0:2^k2-1); %all possible messages with length k2
%add parity bits at the end
dictIII(1:length(dictIII),12)=xor(dictIII(:,1),xor( dictIII(:,2),xor( dictIII(:,4),xor( dictIII(:,5),xor( dictIII(:,7),xor(dictIII(:,9),dictIII(:,11) ))))));
dictIII(1:length(dictIII),13)=xor(dictIII(:,1),xor(dictIII(:,3),xor(dictIII(:,4),xor(dictIII(:,6),xor(dictIII(:,7),xor(dictIII(:,10),dictIII(:,11) ))))));
dictIII(1:length(dictIII),14)=xor(dictIII(:,2),xor(dictIII(:,3),xor(dictIII(:,4),xor(dictIII(:,8),xor(dictIII(:,9),xor(dictIII(:,10),dictIII(:,11) ))))));
dictIII(1:length(dictIII),15)=xor(dictIII(:,5),xor(dictIII(:,6),xor(dictIII(:,7),xor(dictIII(:,8),xor(dictIII(:,9),xor(dictIII(:,10),dictIII(:,11) ))))));

%% |Main Part of Code|-------------------------------------------------------------------------------------------------------------
%Random generate binary messages with equal possibilities per bit
N_msg4=round(randsrc(N,k1,[bits,[0.5,0.5]])); % Nxk array N messages of length k
N_msg11=round(randsrc(N,k2,[bits,[0.5,0.5]]));

for i=1:length(SNR)

  for l=1:N
    %|No code
    no_coding_msg_m=zeros(1,k1);
    no_coding_msg_m=no_coding_msg=N_msg11(l,:);
    no_coding_msg_m(no_coding_msg_m==0)=-1; % modulation 1->"1" and 0->"-1" --bits to symbols

    yI=(sqrt(SNR(i))*no_coding_msg_m)+randn(1,length(no_coding_msg));

    % |HDD|--decision
    yI(yI>0)=1;
    yI(yI<=0)=0;
    %count total ber per message for each snr value(i)
    ber_noCoding(1,i)=ber_noCoding(1,i)+(length(find(yI~=no_coding_msg)))/length(no_coding_msg); %find difference between msg and count

    %% ------------------------------|Hamming (7,4)|--------------------------------------------
    coded_w7=mod(N_msg4(l,:)*G,2); %hamming code by adding parity bits at the end--1xn aka (1x7)
    coded_w7(coded_w7==0)=-1;%bits to symbols ?? modulation
    yII=(sqrt(SNR(i))*coded_w7)+randn(1,n1); %add noise
    % SDD
    for j=1:(k1^2)
      distanceII(j)=norm(yII -dictII(j,:)); % because one vector represent real vector and the other a binary compare its norm
    end
    [minval1,minindx1]=min(distanceII); % find minimum distance and its index in dictionary
    decoded_w7=dictII(minindx1,:); % decode to the word with the minimu difference
    coded_w7(coded_w7==-1)=0; %revert modulation to take back the original msg
    % count ber
    ber_7sdd(1,i)=ber_7sdd(1,i)+(length(find(decoded_w7(1:4)~=coded_w7(1:4))))/n1;

    % HDD
    yII(yII>0)=1;
    yII(yII<=0)=0;
    zII=H*yII';% calculate syndrome
    for g1=1:length(yII)
      if (zII==H(:,g1))
        yII(g1)=~yII(g1);
      end
    end
    ber_7hdd(1,i)=ber_7hdd(1,i)+(length(find(yII~=coded_w7)))/n1;

    %% ------------|Hamming (15,11)|--------------------------------------------------------------
    coded_w15=mod(N_msg11(l,:)*G2,2); % hamming code 1xn --1x15 --adding parity bits to the end
    coded_w15(coded_w15==0)=-1;
    yIII=(sqrt(SNR(i)*coded_w15)+randn(1,n2)); %add noise --chanel output
    % SDD
    for g=1:(k2^2)
      distanceIII(g)=norm(yIII-dictIII(g,:));
    end
    [minval2,minindx2]=min(distanceIII);
    decoded_w15=dictIII(minindx2,:);
    coded_w15(coded_w15==-1)=0; %revert modulation to take original msg to compare it with output
    %count ber
    ber_15sdd(1,i)=ber_15sdd(1,i)+(length(find(decoded_w15(1:11)~=coded_w15(1:11))))/n2;

    % HDD
    yIII(yIII>0)=1;
    yIII(yIII<=0)=0;
    zIII=H2*yIII'; %calculate sydrome
    for g2=1:length(yIII) %15
      if (zIII==H2(:,g2)) %if there is an error fix it
        yIII(g2)=~yIII(g2);
      end
    end
    ber_15hdd(1,i)=ber_15hdd(1,i)+(length(find(yIII~=coded_w15)))/n2;


  end

end

ber_noCoding=ber_noCoding/N; %mean value for all ber
ber_7sdd=ber_7sdd/N;
ber_7hdd=ber_7hdd/N;
ber_15sdd=ber_15sdd/N;
ber_15hdd=ber_15hdd/N;
%% |Plot Results|-----------------------------
figure(1);
semilogy(SNRdB,ber_noCoding);
% set(gca,"XLim", [0 13]);
% set(gca,"YLim", [0 0.5]);
% hold on;
% semilogy(SNRdB,ber_7sdd,"r-s")
% hold on;
% semilogy (SNRdB,ber_7hdd,"k-s")
hold on;
semilogy(SNRdB,ber_15sdd,"k-*")
% hold on;
% semilogy(SNRdB,ber_15hdd,"y-*")
title({"BER for no coding and Hamming","Hamming(7,4),Hamming(15,11)"});
ylabel({"BER"});
xlabel({"SNR(dB)"});
% legend("No coding","Hamming7-SDD","Hamming7-HDD","Hamming15-SDD","Hamming15-HDD")
