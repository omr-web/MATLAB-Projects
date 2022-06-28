clc;
clear all;

% The modulation that will take place according to the mod variable will be selected.
%****************************************************************************
mod='8QAM'; %Select: 8QAM, 16QAM, 64QAM, QPSK, BPSK
%****************************************************************************
switch mod
    case {'BPSK'}
        M=2;
        flag=1;
        Modulation = comm.PSKModulator(2,'BitInput',true,'PhaseOffset',0);
        Demodulation = comm.PSKDemodulator(2,'BitOutput',true,'PhaseOffset',0);
    case {'QPSK'}
        M=4;
        flag=1;
        Modulation = comm.QPSKModulator('BitInput',true);
        Demodulation = comm.QPSKDemodulator('BitOutput',true);
    case {'16QAM'}
        M=16;
        flag=2;
        Modulation = comm.RectangularQAMModulator(16,'BitInput',true);
        Demodulation = comm.RectangularQAMDemodulator(16,'BitOutput',true); 
   case {'64QAM'}
        M=64;
        flag=2;
        Modulation = comm.RectangularQAMModulator(64,'BitInput',true);
        Demodulation = comm.RectangularQAMDemodulator(64,'BitOutput',true);   
   case {'8QAM'}
        M=8;
        flag=2;
        Modulation = comm.RectangularQAMModulator(8,'BitInput',true);
        Demodulation = comm.RectangularQAMDemodulator(8,'BitOutput',true); 
end
%flag değişkeni 
k = log2(M);              % Bits/symbol
numSC = 128;              % Number of OFDM subcarriers (ifft)
cpLen = 32;               % OFDM cyclic prefix length
maxBitErrors = 5e3;       % Maximum number of bit errors
maxNumBits = 1e7;         % Maximum number of bits transmitted

% OFDM modulation  ve demodulation 
ofdmMod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);
ofdmDemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);

%Channel Selection 
%****************************************************************************
Channel_mode = 'RayleighChannel' ; %Select: Rayleigh Channel or 'Rician Channel'
%****************************************************************************
switch Channel_mode
    case {'RicianChannel'}
        Channel = comm.RicianChannel;
        Channel.PathGainsOutputPort =true;
    case {'RayleighChannel'}
         Channel = comm.RayleighChannel;
        Channel.PathGainsOutputPort =true;
end

%White Gauss Noise
awgnChannel = comm.AWGNChannel('NoiseMethod','Variance', 'VarianceSource','Input port');

%Error
errorRate = comm.ErrorRate('ResetInputPort',true);

ofdmDims = info(ofdmMod);

numDC = ofdmDims.DataInputSize(1);
frameSize = [k*numDC 1];%Length of information to be generated from the generator
EbNoVec = (0:10)';%Eb/No
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/numSC); % in 'db' (Eb/No * log2(M)* numDC/numSC)
% Set the SNR vector based on the desired Eb/No range, the number of bits per symbol,... 
% ... and the ratio of the number of data subcarriers to the total number of subcarriers.

berVec = zeros(length(EbNoVec),3);% 1- BER, 2- Total bits transmitted, 3- Number of wrong bits
serVec = zeros(length(EbNoVec),3);% 1- SER, 2- Total bits transmitted, 3- Number of wrong bits
errorStats = zeros(1,3); %Zero is set as the initial value to keep the error values.
for m = 1:length(EbNoVec)
    snr = snrVec(m);
      while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits

          dataIn = randi([0,1],frameSize);                   % Generate binary data
          qpskTx = Modulation(dataIn);                       % Apply  modulation
          txSig = ofdmMod(qpskTx);                           % Apply OFDM modulation
          [fadedSig pathGains] = Channel(txSig);             % Apply fading 
          powerDB = 10*log10(var(fadedSig));                 % Calculate the channel output signal strength
          noiseVar = 10.^(0.1*(powerDB-snr));                % Calculate the noise variance
          rxSig = awgnChannel(fadedSig,noiseVar);            % Pass the signal through a noisy channel
          rxSig_eq=rxSig./sum(pathGains, 2);                 % Equalize (Normalizing the output signal)
          qpskRx = ofdmDemod(rxSig_eq);                      % Apply OFDM demodulation
          dataOut = Demodulation(qpskRx);                    % Apply OFDM demodulation
          errorStats = errorRate(dataIn,dataOut,0);          % Error Calculate
      end
      berVec(m,:) = errorStats;                              % BER data
      serVec(m,:) = errorStats*k;                            % SER data
      errorStats = errorRate(dataIn,dataOut,1);              % Reset the error rate calculator
end

%Drawing theoretical and simulation graphs according to the appropriate modulation technique
if flag==2                                  % For QAM flag=2
    [berTheory,serTheory] = berawgn(EbNoVec,'qam',M);
elseif flag==1                              % For BPSK and QPSK,flag=1
    [berTheory,serTheory] = berawgn(EbNoVec,'psk',M,'nondiff');
end

figure
semilogy(EbNoVec,berVec(:,1),'*')
hold on
semilogy(EbNoVec,berTheory)
title(mod,Channel_mode);

legend('Simulation','Theory','Location','Best')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off

figure
semilogy(EbNoVec,serVec(:,1),'*')
hold on
semilogy(EbNoVec,serTheory)
title(mod,Channel_mode);

legend('Simulation','Theory','Location','Best')
xlabel('Eb/No (dB)')
ylabel('Symbol Error Rate')
grid on
hold off