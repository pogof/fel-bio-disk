format long;

%% Uloha 1
% Komprimuje signál frame-007.bin na bázi kosinové transformace (použijte funkce dct a idct definované v MATLABu). 
% Pro danou kompresi (aproximaci) použijte prvních 75 komponent DCT spektra. Signál je uložen jako binární soubor 
% bez hlavičky, pro načtení do MATLABu použijte funkci loadbin. Původní a dekomprinovaný signál si pro kontrolu ilustrativně zobrazte.

% Spočítejte výkony původního i komprimovaného signálu a určete jaké procento výkonu původního signálu je zahrnuto v signálu komprimovaném.


clear; clc;
x = loadbin('frame-007.bin');
X = dct(x);
 
y = idct(X(1:75), length(X));

Px = mean(x.^2);
Py = mean(y.^2);

pomer = (Py/Px)*100



%% Uloha 2
% Pro posloupnost vzorků [  1  2  3  4  5  6  7  8 ] určete posloupnost doplněnou se sudou 2N symetrií pro účely výpočtu DCT pomocí DFT.
% Správnou odpověď vyberte z následujících možností:

clear; clc;
x = [1 2 3 4 5 6 7 8]; % Zadaná posloupnost

% Určení délky původní posloupnosti
N = length(x);

% Doplnění posloupnosti se sudou 2N symetrií
x_sym = [x fliplr(x)]; % Symetrické doplnění

% Zobrazí původní a doplněnou posloupnost
disp('Původní posloupnost:');
disp(x);
disp('Doplněná posloupnost se sudou 2N symetrií:');
disp(x_sym);



%% Uloha 3
% Spočítejte vyhlazený odhad vzájemné spektrální výkonové hustoty (CPSD) Welchovou metodou pro signály x a y uložené v mat-souboru 
% sig_xy_04.mat (pro načtení do MATLABu použijte "load sig_xy_04.mat"). Signály jsou vzorkované kmitočtem fs = 16 kHz a pro výpočet 
% volte následující parametry:
% 
%     délku krátkodobého segmentu volte 256 vzorků,
%     krátkodobé segmenty váhujte Hammingovým oknem,
%     segmentujte s 50% překryvem,
%     počet bodů FFT volte stejný, jako je délka segmentu,
%     počítejte s implicitním jednostranným odhadem CPSD reálných signálů.
% 
% Určete, který z následujících obrázků je požadovaným odhadem fáze CPSD v radiánech!

clear;clc;
load sig_xy_04.mat;

wlen = 256;
% fs = 16kHz načteno z sig_xy_xx.mat

[Pxy, FF] = cpsd(x, y, hamming(wlen), wlen/2, wlen, fs);

Pxydb = 10*log10(abs(Pxy));

figure;
subplot(211)
plot(FF,20*log10(abs(Pxy)));
title('Magnitude CPSD of selected analyzed frame');
grid
subplot(212)
plot(FF,unwrap(angle(Pxy)))
title('Phase CPSD of selected analyzed frame');
xlabel('f - Frequency [Hz]')
ylabel('phase [rad]')
grid



%% Uloha 4
% Jaký je odstup signálu od šumu (SNR) zašuměného signálu   SX014S01.CS0, je-li referenční čistý signál   SA014S01.CS0? 
% Oba signály jsou uloženy jako binární soubory bez hlavičky, pro načtení do MATLABu použijte funkci loadbin.
% 
% Vypočítané SNR v dB uveďte s přesností na 2 desetinná čísla

clc;clear
noisy = loadbin('SX014S01.CS0');
clean = loadbin('SA014S01.CS0');

noise = noisy - clean;

Pclean = mean(clean.^2);
Pnoise = mean(noise.^2);

SNRdb = 10*log10(Pclean/Pnoise)
disp(['SNRdb: ', num2str(SNRdb)]);



%% Uloha 5
% Určete EUKLIDOVSKOU KEPSTRÁLNÍ VZDÁLENOST na bázi reálného KEPSTRA mezi dvěma signály frame-001.bin a frame-025.bin 
% (oba signály jsou uloženy jako binární soubory bez hlavičky, pro načtení do MATLABu použijte funkci loadbin). Počítejte 
% reálné kepstrum, signály váhujte Hammingovým oknem příslušné délky. Vzdálenost počítejte z prvních 13 koeficientů včetně 
% nultého koeficientu c[0], tj. z koeficientů c[0]-c[12].

% Pro výpočet vzdálenosti použijte funkci cde.m (POZN. Funkci je třeba stáhnout do aktuálního adresáře!!).


% Load the binary signals using loadbin function
s1 = loadbin('frame-001.bin');
s2 = loadbin('frame-025.bin');

% Define parameters
fs = 16000;         % Sampling frequency (assumed)
wlen = 512;         % Frame length
wstep = wlen / 2;   % 50% overlap
cp = 13;            % Number of cepstral coefficients (c[0] to c[12])

% Use a Hamming window for each segment
window = hamming(wlen);

% Compute Real cepstra for both signals
realCepstraRef = vrceps(s1, 1, cp, wlen, wstep, window);
realCepstraDist = vrceps(s2, 1, cp, wlen, wstep, window);

% Select cepstral coefficients c[1]-c[12]
cepstra1 = realCepstraRef(:, 1:13);  % Selecting indices 1 to 13, which corresponds to c[0]-c[12]
cepstra2 = realCepstraDist(:, 1:13);  % Same for the second signal

% Calculate Euclidean cepstral distance for each segment
distances = zeros(size(realCepstraRef, 1), 1);
for i = 1:size(cepstra1, 1)
    distances(i) = cde(cepstra1(i, :), cepstra2(i, :));
end

distance = mean(distances);

% Display the result
disp(['Euclidean cepstral distance: ', num2str(distance)]);



%% Uloha 6
% Určete reálné KEPSTRUM signálu frame-011.bin (uloženo jako binární soubory bez hlavičky, pro načtení do MATLABu 
% použijte funkci loadbin). Signál váhujte Hammingovým oknem příslušné délky. Zobrazte celé vypočítané kepstrum.

x = loadbin('frame-011.bin');
w = hamming(length(x));
x = x.*w;

wlen=512;
channels =1;
cp=16;
wstep=wlen/2;
p=16;   

cr=rceps(x);

figure;
plot(cr);
title('REAL cepstrum of signal')
xlabel('i - frame No. [-]')



%% Uloha 7
% Určete zkreslení delšího signálu SA010S01.CSX na bázi kepstrální vzdálenosti a LPC KEPSTRA, jestliže 
% referenční nezkreslený signál je SA010S01.CS0. Oba signály jsou uloženy jako binární soubory bez hlavičky, 
% pro načtení do MATLABu použijte funkci loadbin. Počítejte LPC kepstrum po segmentech délky wlen=512 s 50% překryvem a 
% uvažujte implicitní váhování každého segmentu Hammingovým oknem. Řád LPC volte p=16, počet kepstrálních koeficentů (bez c[0]) 
% volte cp=20 a vzdálenost počítejte na bázi Euklidovské vzdálenosti včetně nultého koeficientu c[0], tj. z koeficientů  c[0]-c[20].
% Pro výpočet vzdálenosti použijte funkci cde.m (POZN. Funkci je třeba stáhnout do aktuálního adresáře!!).


clear; clc;
refSignal = loadbin('SA010S01.CS0');
distSignal = loadbin('SA010S01.CSX');

% Define parameters
fs = 16000;         % Sampling frequency (assumed)
wlen = 512;         % Frame length
wstep = wlen / 2;   % 50% overlap
p = 16;             % LPC order
cp = 20;            % Number of cepstral coefficients (c[0] to c[20])

% Use a Hamming window for each segment
window = hamming(wlen);

% Compute LPC cepstrum for both signals
lpcCepstraRef = vaceps(refSignal, 1, p, cp, wlen, wstep, window);
lpcCepstraDist = vaceps(distSignal, 1, p, cp, wlen, wstep, window);

% Compute Euclidean cepstral distance for each frame and average
numFrames = size(lpcCepstraRef, 1);
distances = zeros(numFrames, 1);
for i = 1:numFrames
    c1 = lpcCepstraRef(i, :);
    c2 = lpcCepstraDist(i, :);
    distances(i) = cde(c1, c2);
end

% Average distance
avgDistance = mean(distances);

% Display the result
disp(['Average Euclidean cepstral distance: ', num2str(avgDistance)]);



%% Uloha 8
% Pro signály sig1 a sig2 vzorkované kmitočtem fs = 16 kHz a uložené v mat-souboru sigs_2chan_02.mat (pro načtení do MATLABu použijte "load sigs_2chan_02.mat") vypočtěte koherenční funkci, konkrétně MSC (Magnitude Square Coherence), přičemž pro výpočet volte následující parametry:
% 
%     délka krátkodobého segmentu - 16 ms,
%     váhování - Hammingovo okno odpovídající délky,
%     segmentace - s 75% překryvem,
%     řád FFT - stejný, jako je délka krátkodobého segmentu.
% 
% Určete průměrnou koherenci (tj. průměrnou hodnotu vypočítané MSC). Výsledek uveďte s minimální přesností na 3 platné cifry.

clear; clc;
load sigs_2chan_02.mat;

fs = 16000; % Kmitočet vzorkování
wlen = 0.016 * fs; % Délka krátkodobého segmentu v počtu vzorků
overlap = 0.75; % Překryv segmentů (75%)
noverlap = overlap * wlen; % Počet vzorků překryvu
nfft = wlen; % Řád FFT

% Vypočítání koherenční funkce
[Cxy, F] = mscohere(sig1, sig2, hamming(wlen), noverlap, nfft, fs);

% Výpočet průměrné koherence
mean_Cxy = mean(Cxy);

disp(['Průměrná koherence: ', num2str(mean_Cxy)]);







