%% TSKS10 Laborationskod
% Gustav Jannering
% gusja113

% Återställ
close all;
clear;

% Läs in ljudfilen
[signal, sampel_frek] = audioread('signal-gusja113.wav');
signal_l = length(signal);

% Transformera signalen, 
fft_signal = fft(signal);

% Halvera signalen för att undvika vikning
fft_signal = fft_signal(1:signal_l/2);

% Beräkna vektorer som kommer att användas som axlar i grafer
f_axel = (0:signal_l/2-1)*sampel_frek/signal_l;
t_axel = linspace(0,signal_l/sampel_frek, signal_l);

% Plotta |Y(f)| för att undersöka signalen i frekvsensspektrumet
figure;
plot(f_axel, abs(fft_signal));
title('Amplitud spektrum');
xlabel('Frekvens f [kHz]');
ylabel('Amplitud');


% Bandbredden är 20 kHz
bandbredd = 20000;

% Möjliga värden på fc, från |Y(f)|
frek_1 = 36000; % 36 kHz
frek_2 = 74000; % 74 kHz
frek_3 = 112000; % 112 kHz

%skapa ny figur
figure;

% Bandpassfilter för att undersöka första möjliga fc, 36 kHz 
[B1,A1] = butter(10, [30000,50000]/(sampel_frek/2));

delsignal1 = filter(B1, A1, signal);

subplot(3,1,1);
plot(t_axel, delsignal1);
xlabel('Tid [s]');
ylabel('Amplitud');
title('Bärfrekvens 36 kHz');

% Bandpassfilter för att undersöka andra möjliga fc, 74 kHz 
[B2,A2] = butter(10, [60000,80000]/(sampel_frek/2));

delsignal2 = filter(B2, A2, signal);

subplot(3,1,2);
plot(t_axel, delsignal2);
xlabel('Tid [s]');
ylabel('Amplitud');
title('Bärfrekvens 74 kHz');

% Bandpassfilter för att undersöka tredje möjliga fc, 112 kHz 
[B3,A3] = butter(10, [110000,120000]/(sampel_frek/2));

delsignal3 = filter(B3, A3, signal);

subplot(3,1,3);
plot(t_axel, delsignal3);
xlabel('Tid [s]');
ylabel('Amplitud');
title('Bärfrekvens 112 kHz');

% Bandpassfilter som används för att bestämma f1,f2
[B4,A4] = butter(10, [100000,103000]/(sampel_frek/2));

delsignal4 = filter(B4, A4, signal); %delsignal som innehåller w(t) 
del_l = length(delsignal4); % w(t)s längd
fft_del4 = fft(delsignal4);

% Halvera signalen för att undvika vikning
fft_del4 = fft_del4(1:del_l/2);

f_axel_ = sampel_frek/2*linspace(0, 1, del_l/2);

figure;
plot(f_axel, abs(fft_del4));

fc = 112000; % Från |Y(t)| 
f1 = 102500; % Från |Y(t)|
f2 = 102501; % Från |Y(t)|

% Autokorrelation för att undersöka tidsfördröjningen
[korrelation,lags] = xcorr(delsignal2);

%Endast positiva värden är intressanta
lags = lags(lags>0);
korrelation = korrelation(lags>0);

figure;
plot(lags/sampel_frek, abs(korrelation));
xlabel('Tid [s]');
title('Korrelation');

% Från figur
fordojning = 0.39;
sampel_fordroj = fordojning * sampel_frek;
sampel_f_l = size(sampel_fordroj);

eko_signal = zeros(sampel_f_l);

% Punkter före eko uppstår
for n = 1 : sampel_fordroj
   
    eko_signal(n) = delsignal3(n);

end
% Punkter efter eko uppstår
% Ta bort eko
for n = sampel_fordroj+1 : signal_l
   
    eko_signal(n) = delsignal3(n) - 0.9*eko_signal(n-sampel_fordroj);

end

% LP-filter
[B5,A5] = butter(10,bandbredd/(sampel_frek/2),'low');

% Fasvridningskonstant
phi = 0.8; %pi/4 < phi < pi/3 

for n = 0: 0.1:pi/2
    disp(n)
    test = 2*pi*fc*t_axel.'+n;

    %I/Q-demodulering
    X_i = filter(B5, A5, eko_signal.*(2*cos(test.')));
    
    X_q = -filter(B5, A5, eko_signal.*(2*sin(test.')));
    
    %Downsampla med 10 för att kunna lyssna
    i = downsample(X_i, 10);
    q = downsample(X_q, 10);
    
    %Uppspelning
    %soundsc(i, sampel_frek/10);
    %soundsc(q, sampel_frek/10);
    pause;
end


