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
        first_part=117;
        last_part=234;
        Modulation = comm.PSKModulator(2,'BitInput',true,'PhaseOffset',0);
        Demodulation = comm.PSKDemodulator(2,'BitOutput',true,'PhaseOffset',0);
    case {'QPSK'}
        M=4;
        flag=1;
        first_part=234;
        last_part=468;
        Modulation = comm.QPSKModulator('BitInput',true);
        Demodulation = comm.QPSKDemodulator('BitOutput',true);
    case {'16QAM'}
        M=16;
        flag=2;
        first_part=468;
        last_part=936;
        Modulation = comm.RectangularQAMModulator(16,'BitInput',true,....
    'NormalizationMethod','Average power','SymbolMapping','Gray');
        Demodulation = comm.RectangularQAMDemodulator(M,'BitOutput',true,....
    'NormalizationMethod','Average power','SymbolMapping','Gray'); 
   case {'64QAM'}
        M=64;
        flag=2;
        first_part=702;
        last_part=1404;
        Modulation = comm.RectangularQAMModulator(64,'BitInput',true,....
    'NormalizationMethod','Average power','SymbolMapping','Gray');
        Demodulation = comm.RectangularQAMDemodulator(64,'BitOutput',true,....
    'NormalizationMethod','Average power','SymbolMapping','Gray');  
    
   case {'8QAM'}
        M=8;
        flag=2;
        first_part=351;
        last_part=702;
        Modulation = comm.RectangularQAMModulator(8,'BitInput',true,....
    'NormalizationMethod','Average power','SymbolMapping','Gray');
        Demodulation = comm.RectangularQAMDemodulator(8,'BitOutput',true,....
    'NormalizationMethod','Average power','SymbolMapping','Gray'); 
end

k = log2(M);           % Bits/symbol
numSC = 128;           % Number of OFDM subcarriers (ifft)
cpLen = 32;            % OFDM cyclic prefix length
maxBitErrors = 5e3;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted


ofdmMod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen); %OFDM Modulator
ofdmDemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen); %OFDM Demodulator

%Channel Selection 
%****************************************************************************
Channel_mode = 'Rician' ; %Select: 'Rayleigh' or 'Rician'
%****************************************************************************

mimochannel = comm.MIMOChannel( ...
    'TransmitCorrelationMatrix' ,[1 1;1 1], ... 
    'ReceiveCorrelationMatrix' ,[1 1;1 1], ...
    'NumTransmitAntennas',2,'NumReceiveAntennas',2, ...
    'FadingDistribution',Channel_mode, ...
    'PathGainsOutputPort', true);                         % 2x2 MIMO Channel


awgnChannel = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');                       %AWGN Noise


errorRate = comm.ErrorRate('ResetInputPort',true);    %Error Rate Calculation

ofdmDims = info(ofdmMod); %Use the 'info' object function to get the modulator input data

numDC = ofdmDims.DataInputSize(1);
frameSize = [k*numDC 1]; %Length of information to be generated from the generator


EbNoVec = (0:10)'; %Eb/No
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/numSC); % in 'db' (Eb/No * log2(M)* numDC/numSC)
% Set the SNR vector based on the desired Eb/No range, the number of bits per symbol,... 
% ... and the ratio of the number of data subcarriers to the total number of subcarriers.

berVec = zeros(length(EbNoVec),3); % 1- BER, 2- Total bits transmitted, 3- Number of wrong bits
serVec = zeros(length(EbNoVec),3);% 1- SER, 2- Total bits transmitted, 3- Number of wrong bits
errorStats = zeros(1,3); %Zero is set as the initial value to keep the error values.

for m = 1:length(EbNoVec)
    snr = snrVec(m);
    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
        dataIn1 = randi([0,1],frameSize);                                % Generate binary data
        dataIn = [dataIn1 ; dataIn1];                                    % Data multiplexing
        qamTx1 = Modulation(dataIn(1:first_part));                       % Apply QAM modulation
        qamTx2 = Modulation(dataIn(first_part+1:last_part));             % Apply QAM modulation
        txSig1 = ofdmMod(qamTx1);                                        % Apply OFDM modulation
        txSig2 = ofdmMod(qamTx2);                                        % Apply OFDM modulation
        txSig = [txSig1 txSig2];                                         % Concatenate

        [fadedSig pathGains] = mimochannel(txSig);    % Apply fading and get the path gains
        
        pathGains1_1=pathGains(:,:,1,1);              % Path-1
        pathGains2_1=pathGains(:,:,2,1);              % Path-2
        pathGains1_2=pathGains(:,:,1,2);              % Path-3
        pathGains2_2=pathGains(:,:,2,2);              % Path-4

        powerDB = 10*log10(var(fadedSig));            % Calculate the channel output signal strength
        noiseVar = 10.^(0.1*(powerDB-snr));           % Calculate the noise variance
        rxSig = awgnChannel(fadedSig,noiseVar);       % Pass the signal through a noisy channel

        pathGains1=[pathGains1_1 pathGains2_1];       %1. Transmitter paths
        pathGains2=[pathGains1_2 pathGains2_2];       %2. Transmitter paths
        rxSig_eq1=rxSig./sum(pathGains2, 2);          % Equalize (Normalizing the output signal)
        rxSig_eq2=rxSig./sum(pathGains1, 2);          % Equalize (Normalizing the output signal)
          
        path1 = rxSig_eq1(:,1);                       %Path selection
        path2 = rxSig_eq2(:,2);                       %Path selection
        qamRx1 = ofdmDemod(path1);                    % Apply OFDM demodulation
        qamRx2 = ofdmDemod(path2);                    % Apply OFDM demodulation
        
        dataOut1 = Demodulation(qamRx1);              % Apply QAM demodulation
        dataOut2 = Demodulation(qamRx2);              % Apply QAM demodulation
        dataOut = [dataOut1 ; dataOut2];
        errorStats = errorRate(dataIn,dataOut,0);     % Collect error statistics
    end
    berVec(m,:) = errorStats;                              % BER data
    serVec(m,:) = errorStats*k;                            % SER data
    errorStats = errorRate(dataIn,dataOut,1);         % Reset the error rate calculator
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