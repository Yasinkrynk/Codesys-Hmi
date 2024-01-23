% MIDTERM-HW13 CMS SOFT NONLINEAR
% Doç. Dr. Sela_mi Beyhan
clc, clear all, close all

%% Zaman değerleri---PLC de 5sn yeterli
T  = 20;       % Simülasyon süresi
Ts = 0.001;    % Örnekleme adımı

%% Durumlar ve girişlerin ilk degerleri
x = zeros(2,200);
u = zeros(1,200);

%% Durumların başlangıç değerleri
x(:,1) = [0;0];    % Başlangıç durumu
xref(:,1) = [2;2]; % Referans başlangıç durumu

%% System Matrices
a11 = 0; a12 = 1; a21 = -3; a22 = -2;
b21 = 3; c11 = 0; c12 = 2;  
%b11 = 0; SMC için sistem mecbur oyle x1dot = x2 seklınde

%% PID Kontrol sinyali
error = [0;0];             % Bunun daha sonra indekslenmesine gerek yok
umax  = 12;                % Her iki u degeri sınırlanacak
ks    = 0.01;              % US sabiti
lamda = 50;                % SS sabiti
Thres = 0.1;               % Saturasyon sabiti

%% Ana Döngü: PID kontrol ile kontrol işareti üret ve sistemi Euler ile integre et
for k = 1:(T/Ts)
    t(k) = k*Ts;                      % Sürekli zaman
    
    cc = 1; ww = 2; aa = 2.5;
    xrefdot(:,k) = [xref(2,k); cc*(1-xref(1,k)^2)*xref(2,k)-xref(1,k)+aa*cos(ww*t(k))];  
    xref(:,k+1) = xref(:,k) + Ts*xrefdot(:,k);   

 
    error(1,k) = xref(1,k)-x(1,k);     % Referans 1 hatası
    error(2,k) = xref(2,k)-x(2,k);     % Referans 2 hatası
    surface(k) = error(2,k)+error(1,k)*lamda;
    
    if abs(surface(k))<Thres
        sat(k) =  ks*error(1,k);
    else
        sat(k) = -ks*tanh(error(1,k));
    end
    
    % Kontrol sinyali türetme
    fhat = a21*x(1,k)+a22*x(2,k);  
    ghat = b21;
    if k>2
        dot_deriv(k) = (xref(2,k)-xref(2,k-1))/Ts;
    else 
        dot_deriv(k) = 0;
    end
    
    u(k) = (-fhat+surface(k)+dot_deriv(k)+sat(k))/ghat;
    
    if (u(k)>umax),  u(k) =  umax; end
    if (u(k)<-umax), u(k) = -umax; end
    
    % PLC de x1 ve x2 integrasyon ayrı ayrı devam ediyor burada vektörel verilde
    xdot(:,k) = [a11*x(1,k)+a12*x(2,k);a21*x(1,k)+a22*x(2,k)+b21*u(k)];  % Vektörel xdot
    x(:,k+1) = x(:,k) + Ts*xdot(:,k);                                    % Vektörel x integrasyonu
end
%% PLC ile karşılaştırılacak veriler
x1_verisi = x(1,1:20)
x2_verisi = x(2,1:20)
u_verisi  = u(1:20)

%% Grafikleri çiz
close all
figure;
plot(t, xref(1,1:end-1), '-k', 'linewidth', 1.5); % Referans sinyali
hold on;
plot(t, x(1,1:end-1), '--b', 'linewidth', 2); % Sistem çıkışı
grid on;
legend('Referans', 'THI-1','Location','southeast');
legend('boxoff')
xlabel('Zaman (s)');
ylabel('x1: 1. Sistem Çıkışı');

figure;
plot(t, xref(2,1:end-1), '-k', 'linewidth', 1.5); % Referans sinyali
hold on;
plot(t, x(2,1:end-1), '--b', 'linewidth', 2); % Sistem çıkışı
grid on;
legend('Referans', 'THI-2','Location','southeast');
legend('boxoff')
xlabel('Zaman (s)');
ylabel('x2: 2. Sistem Çıkışı');

figure;
plot(t, u(1,:), '-r', 'linewidth', 1.5); % Kontrol sinyali
grid on;
xlabel('Zaman (s)');
ylabel('Kontrol sinyali');



