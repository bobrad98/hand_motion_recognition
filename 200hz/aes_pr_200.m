clear
close all
clc

%% unos podataka, definisanje konstanti

%obucavajuci
data1 = load('ob1.txt');       %dodir prstima
data2 = load('ob2.txt');       %pesnica
data3 = load('ob3.txt');       %spolja
data4 = load('ob4.txt');       %unutra

%testirajuci
data11 = load('t1.txt');
data22 = load('t2.txt');
data33 = load('t3.txt');
data44 = load('t4.txt');

Ts = 0.005;               %vremenski razmak izmedju odbiraka podataka
fs = 1/Ts;               %frekvencija akvizicije
cut = 0.5;

N = 140;     %sirina prozora
del = Ts;   %razmak izmedju susednih prozora

%% filtriranje, dodavanje informacije o klasi

data1 = fil(data1, fs, cut);
data2 = fil(data2, fs, cut);
data3 = fil(data3, fs, cut);
data4 = fil(data4, fs, cut);

Xo1 = obelezja(data1, N);
Xo2 = obelezja(data2, N);
Xo3 = obelezja(data3, N);
Xo4 = obelezja(data4, N);

Xo1(:, end + 1) = 1;
Xo2(:, end + 1) = 2;
Xo3(:, end + 1) = 3;
Xo4(:, end + 1) = 4;

Xo = [Xo1; Xo2; Xo3; Xo4];
[odbo, atr] = size(Xo);

%%%%%%

data11 = fil(data11, fs, cut);
data22 = fil(data22, fs, cut);
data33 = fil(data33, fs, cut);
data44 = fil(data44, fs, cut);

Xt1 = obelezja(data11, N);
Xt2 = obelezja(data22, N);
Xt3 = obelezja(data33, N);
Xt4 = obelezja(data44, N);

Xt1(:, end + 1) = 1;
Xt2(:, end + 1) = 2;
Xt3(:, end + 1) = 3;
Xt4(:, end + 1) = 4;

Xt = [Xt1; Xt2; Xt3; Xt4];
[odbt, ~] = size(Xt);

%% neuralna

%obucavanje
Xn = Xo(:, 1:end-1)';      %ulazni vektori
Tn = zeros(4, odbo);
Tn(1, 1:odbo/4) = 1;
Tn(2, odbo/4+1:odbo/2) = 1;
Tn(3, odbo/2+1:3*odbo/4) = 1;
Tn(4, 3*odbo/4 + 1: end) = 1;      %zeljeni izlaz

%net = newff(Xn, Tn, (atr - 1)/2, {'tansig', 'logsig'});
net = patternnet((atr-1)/2);
net.divideFcn = '';
%{
net.divideParam.trainRatio = 0.85;  
%net.divideParam.valRatio = 0.15;
net.divideParam.testRatio = 0.15;
%}

%testiranje
Xnt = Xt(:, 1:end-1)';
Tnt = zeros(4, odbt);
Tnt(1, 1:odbt/4) = 1;
Tnt(2, odbt/4+1:odbt/2) = 1;
Tnt(3, odbt/2+1:3*odbt/4) = 1;
Tnt(4, 3*odbt/4 + 1: end) = 1;      %zeljeni izlaz

y = train(net, Xn, Tn);
P = y(Xnt);
%p = zeros(length(P));
%{
for k = 1:length(P)
    for j = 1:4
        if
    end
end
%}
%C = confusionmat(Xnt(:, end), P);
%cmat(i) = C(end, end);
plotconfusion(Tnt, P);
%P = y(Xnt);
%P = round(P);
%plotconfusion(Tnt, P);

%pusti u for loop 20 komada pa onda nadji stdev

function dataf = fil(data, fs, cut)
    [l, sig] = size(data);
    dataf = zeros(l, sig);
    %ts = 1/fs;
    %t = 0:l:l*ts;    %vremenska osa
    [b, a] = butter(4, 2*cut/fs);
    for i = 1:sig
        dataf(:, i) = filter(b, a, data(:, i));
    end
end

function X = obelezja(data, N)
    dim = max(size(data));
    Wn = dim - N;               %broj prozora
    X = zeros(Wn, min(size(data)) * 3);
    for i = 1:8
        for j = 1:Wn
            data1 = data(j : j + N - 1, i);     %podaci u trenutnom prozoru
            data2 = zeros(1, length(data1) - 1);
            for k = 1:length(data1) - 1
                data2(k) = data1(k + 1)  - data1(k);
            end
            X(j, 3*(i - 1) + 1) = sum(abs(data1)) / N;          %mav
            X(j, 3*(i - 1) + 2) = sqrt(sum(data1 .* data1) / N); %rms
            X(j, 3*i) = sum(abs(data2));  %WL
        end   
    end
end
