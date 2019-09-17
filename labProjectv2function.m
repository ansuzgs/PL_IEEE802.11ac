function [numErrors] = labProjectv2function(AWGN,snr, c, yQPSKofdm_CP, Nifft, Nsym, CPSize, tau, ts, Ndata, hMQAM, x, k)

    %% Canal
    reset(c);
    [yQPSKrx, pathGains] = c(yQPSKofdm_CP);
    yQPSK_AWGN = step(AWGN, yQPSKrx);
    
    %% Demodulación OFDM
    yrxconCP = reshape(yQPSK_AWGN, Nifft+CPSize, Nsym);
    yrxsinCP = zeros(Nifft, Nsym);
    yrxsinCP = yrxconCP(CPSize+1:end, :);

    zrx = fft(yrxsinCP, Nifft);
    
    %% Ecualización
    h = zeros(max(tau/ts)+1,1);
    h(1) = pathGains(1,1);
    h(6) = pathGains(1,2);
    h(8) = pathGains(1,3);
    h(11) = pathGains(1,4);
    Hest = fft(h,Nifft);
    zrx = zrx./Hest;
    
    %% Obtención trama M-QAM
    tramaMQAMrx = zeros(Ndata,Nsym);
    tramaMQAMrx(1:9,:) = zrx(2:10,:);
    tramaMQAMrx(10:18,:) = zrx(12:20,:);
    tramaMQAMrx(19:26,:) = zrx(22:29,:);
    tramaMQAMrx(27:34,:) = zrx(37:44,:);
    tramaMQAMrx(35:43,:) = zrx(46:54,:);
    tramaMQAMrx(44:end,:) = zrx(56:end,:);
    
    %% Demodulación M-QAM
    ydemo = step(hMQAM, tramaMQAMrx(:));
    
    %% Cálculo de errores
    numErrors = sum(abs(x-ydemo));
    
end