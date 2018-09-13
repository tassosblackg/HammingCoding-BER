%Final project theoria plhroforias--Hamming Coding,BER
%Implementation on Octave 4.2.2 on Ubuntu 16.04

%  |Warning: if you run it at Matlab you have to change xor()

% *|NOTICE: parity bits are added at the end of word , not at their regular position
% *|e.g word p1p2p3..(1001 001), word=1001 parity=001 p1=0,p2=0,p3=1 --from left to right
% *|suppose adding noise will cause error always

clear;
clc;

% need this Octave
pkg load communications;
pkg load signal;
%

bits=[1 0];
N=10000; %one million recursions montecarlo
SNR=[0:3:12]; %SNR vector repeat experiment for different SNR
%SNR=10.^(SNRdB./10)
% init ber vectors
ber_snr7i=zeros(1,length(SNR));
ber_snr15i=zeros(1,length(SNR));
ber_snr7ii=zeros(1,length(SNR));
ber_snr15iii=zeros(1,length(SNR));

% ------Hamming(7,4)
%G add parity bits in the 2^i position
%G=[1 1 0 1; 1 0 1 1; 1 0 0 0; 0 1 1 1; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %generator matrix
G=[1 0 0 0 0 1 1; 0 1 0 0 1 0 1; 0 0 1 0 1 1 0; 0 0 0 1 1 1 1]; %generator matrix ,coded_word=word*G--parity bits at the end
H=[0 1 1 1 1 0 0; 1 0 1 1 0 1 0; 1 1 0 1 0 0 1];
%H=[1 0 1 0 1 0 1; 0 1 1 0 0 1 1; 0 0 0 1 1 1 1]; %z=H*r sydrome
%------Hamming(15,11)
%G2 add parity bits at the end of word
G2=[1 0 0 0 0 0 0 0 0 0 0 1 1 0 0; 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0;
    0 0 0 1 0 0 0 0 0 0 0 1 1 1 0; 0 0 0 0 1 0 0 0 0 0 0 1 0 0 1; 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1
    0 0 0 0 0 0 1 0 0 0 0 1 1 0 1; 0 0 0 0 0 0 0 1 0 0 0 0 0 1 1; 0 0 0 0 0 0 0 0 1 0 0 1 0 1 1
    0 0 0 0 0 0 0 0 0 1 0 0 1 1 1; 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1];
H2=[1 0 1 1 1 0 0 0 1 1 1 1 0 0 0; 1 1 0 1 1 0 1 1 0 0 1 0 1 0 0;
    1 1 1 0 1 1 0 1 1 0 0 0 0 1 0; 1 1 1 1 0 1 1 0 0 1 0 0 0 0 1];

for l=1:length(SNR)
  tBER7=tBER15=tBERii=tBERiii=0; %total BER from all recursions #N
  for j=1:N
    %init stage) information bits
    %for I)
    msg1=round(randsrc(1,4,[bits,[0.5,0.5]]));
    msg2=round(randsrc(1,11,[bits,[0.5,0.5]]));


    % II)-hamming(7,4)
    coded_w7=mod(msg1*G,2); %create coded word 1x7

    % III)-hamming(15,11)
    coded_w15=mod(msg2*G2,2); %create coded word 1x15

    %create white gaussian noise [0-12]dB with step 3dB
    % n1=awgn(1,7,SNR(l),'dB');
    % n2=awgn(1,15,SNR(l),'dB');

    %apply white gaussian noise to msg through channel--NO CODING
    ch_out1=awgn(msg1,SNR(l));
    ch_out2=awgn(msg2,SNR(l));

    %---hamming add noise
    ch_outII=awgn(coded_w7,SNR(l)); %chanel output of msg with hamming code (7,4)
    ch_outIII=awgn(coded_w15,SNR(l));

    %---------|NO CODING|---------------------------------------
    %% --Ι)
    %Output of NOT encode message
    y7i=zeros(1,length(ch_out1));
    y15i=zeros(1,length(ch_out2));

    %decision for 7bits
    y7i(ch_out1>0)=1;
    y7i(ch_out1<=0)=0;
    %decision for 15bits
    y15i(ch_out2>0)=1;
    y15i(ch_out2<=0)=0;

    %calculate BER i) for each experement--one loop
    BER7=(sum(abs(ch_out1-y7i)))/length(ch_out1);
    BER15=(sum(abs(ch_out2-y15i)))/length(ch_out2);

    tBER7=tBER7+BER7; %total BER for N
    tBER15=tBER15+BER15;
    %------------------------------------------------------------------

    % HAMMING Coding --decide from signal with noise
    yII(ch_outII>0)=1;
    yII(ch_outII<=0)=0;

    yIII(ch_outIII>0)=1;
    yIII(ch_outIII<=0)=0;

    % calculate sydrome
    %if error-> position-column of H where error
    zII=H*yII';
    zIII=H2*yIII';

    %correct the error
    %--II)
    for g=1:7
      if(zII==H(:,g)) %
        yII(g)=~yII(g); %flip the error bit  to fix
      end
    end
    %III)
    for g=1:11
      if(zIII==H2(:,g)) %
        yIII(g)=~yIII(g); %flip the error bit to fix
      end
    end

    %CREATE DICTIONARY FOR SDD
    %II)
    dictII=de2bi(0:2^4-1);
    %add parity bits at the end
    dictII(1:length(dictII),5)=xor(dictII(:,1),xor(dictII(:,2),dictII(:,4)));
    dictII(1:length(dictII),6)=xor(dictII(:,1),xor(dictII(:,3),dictII(:,4)));
    dictII(1:length(dictII),7)=xor(dictII(:,2),xor(dictII(:,3),dictII(:,4)));

    %III)
    dictIII=de2bi(0:2^11-1);
    %add parity bits at the end
    dictIII(1:length(dictIII),12)=xor(dictIII(:,1),dictIII(:,2),dictIII(:,4),dictIII(:,5),dictIII(:,7),dictIII(:,9),dictIII(:,11));
    dictIII(1:length(dictIII),13)=xor(dictIII(:,1),dictIII(:,3),dictIII(:,4),dictIII(:,6),dictIII(:,7),dictIII(:,10),dictIII(:,11));
    dictIII(1:length(dictIII),14)=xor(dictIII(:,2),dictIII(:,3),dictIII(:,4),dictIII(:,8),dictIII(:,9),dictIII(:,10),dictIII(:,11));
    dictIII(1:length(dictIII),15)=xor(dictIII(:,5),dictIII(:,6),dictIII(:,7),dictIII(:,8),dictIII(:,9),dictIII(:,10),dictIII(:,11));

    % --| Minimum distance|aka SDD
    % II)
    for m=1:(2^4)           %calculate distance of received word with every possible word--dictII
      distanceII(m)=norm(yII-dictII(m,:));
    end
    [minval,minindx]=min(distanceII);       %Finding index of the minimum distance
    decoded_wII=dictII(minindx,:); %choose word with minimum distance @@

    tBERii=tBERii+(sum(abs(coded_w7-decoded_wII)))/length(decoded_wII); %calculate error per coded word, and sum all of them

    % III)
    for m=1:(2^11)           %calculate distance of received word with every possible word--dictII
      distanceIII(m)=norm(yIII-dictIII(m,:));
    end
    [minval2,minindx2]=min(distanceIII);       %Finding index of the minimum distance
    decoded_wIII=dictIII(minindx2,:); % @@

    tBERiii=tBERiii+(sum(abs(coded_w15-decoded_wIII)))/length(coded_w15);

  end


    % BER case i)
    ber_snr7i(l)=tBER7/N;
    ber_snr15i(l)=tBER15/N;

    % BER case ii)
    ber_snr7ii(l)=tBERii/N;

    % BER case iii)
    ber_snr15iii(l)=tBERiii/N;
end
%-----------------------------------------------

figure(1);
semilogy(SNR,ber_snr7i,"b-<")
hold on;
semilogy(SNR,ber_snr15i,"b->")
hold on;
plot(SNR,ber_snr7ii,"r-*")
hold on;
plot(SNR,ber_snr15iii,"g-s")
title({"BER for no coding and Hamming","Hamming(7,4),Hamming(15,11)"});
ylabel({"BER"});
xlabel({"SNR(dB)"});
legend("No Coding 7","No coding 15","Hamming 7","Hamming 15")
%text('')
