# HammingCoding-BER
Project for testing BER on hamming code-Information Theory subject

This is a matlab script writen in octave on linux


The purpose of this src code is to estimate BER (Bit Error Rate) while changing SNR.

We have a random generated binary message that is send through a channel where WhiteGaussian Noise is applied

And we compare the received message to the original while counting the bit error between two sequences.

Comparasion is made through three cases :
-------------------------------------------

-1st Message without any coding 
  -Message length= 4
  -Message length= 11
-2nd Message is coded with Hamming(7,4)
-3d  Message is code with Hamming(15,11)

Decoding :
--------------

In case 1, we use HDD(Hard Decision Decoding) but in case 2 & 3 we use SDD(Soft Decision Decoding)

Finally plot the BER(SNR) function for each case.
