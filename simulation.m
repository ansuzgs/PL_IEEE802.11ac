%% Parámetros globlales
m = [4,16,64]; % modulaciones
snr = 10:2:30; % Energy per bit to noise ratio (dB)
Nit = 200; % Número de iteraciones de simulación
Nbits = 1e5; % Tamaño inicial de datos
Ndata = 52; % Portadoras de datos
Nifft = 64; % Tamaño de la FFT
Nnull = 8; % Portadoras nulas
Npiloto = 4; % Portadoras piloto
fs = 20e6; % Ancho de banda
tcp = 800e-9; % Duración banda de guarda
ts = 1/fs; % Tiempo de simbolo
fd = 1; % Frecuencia doppler
tau = [0, 5, 7, 10].*ts; % Retardos del canal multicamino
G = [0, -3, -5, -10]; % Ganancias del canal multicamino


%% Simulación
for iii = 1:3
    switch(iii)
        case 1
            disp('QPSK')
        case 2
            disp('16-QAM')
        case 3
            disp('64-QAM')
        otherwise
            disp('Error')
    end
    
    ber = zeros(length(snr), 1); % Inicialización BER
    errores2 = zeros(Nit, 1); % Inicialización Errores
    M = m(iii); % Orden de la modulación [4, 16, 64]
    k = log2(M); % Bits por simbolo
    Nsym = ceil(Nbits/(log2(M)*Ndata)); % Número de simbolos necesarios
    Nbits2 = Nsym*k*Ndata; % Numero de bits necesarios


    %% Generacion de los bits
    x = randi([0 1], Nbits2, 1);

    %% M-QAM
    hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true);
    dataMod = step(hMod,x);
    datos = reshape(dataMod, Ndata, Nsym);
    piloto = step(hMod, ones(k,1));

    %% Preparación trama OFDM
    tramaOFDM = zeros(Nifft,Nsym);
    tramaOFDM(2:10,:) = datos(1:9,:);
    tramaOFDM(11,:) = piloto;
    tramaOFDM(12:20,:) = datos(10:18,:);
    tramaOFDM(21,:) = piloto;
    tramaOFDM(22:29,:) = datos(19:26,:);
    tramaOFDM(37:44,:) = datos(27:34,:);
    tramaOFDM(45,:) = piloto;
    tramaOFDM(46:54,:) = datos(35:43,:);
    tramaOFDM(55,:) = piloto;
    tramaOFDM(56:end,:) = datos(44:end,:);

    %% Modulación OFDM
    yQPSKofdm = ifft(tramaOFDM,Nifft); 

    %% Introducir CP
    CPSize = tcp/ts;
    yQPSKofdm_CP = zeros(Nifft+CPSize,Nsym);
    yQPSKofdm_CP(1:CPSize,:) = yQPSKofdm(end-(CPSize-1):end,:);
    yQPSKofdm_CP(CPSize+1:end,:) = yQPSKofdm(:,:);

    %% Canal
    yQPSKofdm_CP = yQPSKofdm_CP(:);
    c = comm.RayleighChannel('SampleRate', fs, 'PathDelays', tau, 'AveragePathGains', G, 'MaximumDopplerShift', fd, 'PathGainsOutputPort',true);
    AWGN = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
    pot = mean(abs(yQPSKofdm_CP).^2);
    AWGN.SignalPower = pot;

    %% Definición demodulador
    hMQAM = comm.RectangularQAMDemodulator('ModulationOrder',M, 'BitOutput', true);


    %% Simulación para cada modulació
    for i =1:length(snr)
        disp(['SNR = ', num2str(snr(i)), 'dB'])
        errores2 = zeros(Nit, 1);
        AWGN.SNR = snr(i);
        for ii = 1:Nit
            [errores2(ii)] = labProjectv2function(AWGN, snr(i), c, yQPSKofdm_CP, Nifft, Nsym, CPSize, tau, ts, Ndata, hMQAM, x, k);
        end
        ber(i) = sum(errores2)/(Nit*Nbits2);
    end

    semilogy(snr, ber, '-o')
    grid on
    xlabel('SNR')
    ylabel('BER')
    hold on
end

legend('QPSK','16-QAM','64-QAM')
title('Comparación BER-SNR')
