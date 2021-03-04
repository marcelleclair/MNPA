set(0, 'DefaultFigureWindowStyle', 'docked');
close all

% Component Values
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000;
C1 = 0.25;
L = 0.2;
alpha = 100;

C = zeros(6,6);
C(2,1) = C1;
C(2,2) = -C1;
C(6,6) = L;

G = zeros(6,6);
G(1,1) = 1;
G(2,1) = 1/R1;
G(2,2) = -(R1+R2)/(R1*R2);
G(2,6) = -1;
G(3,3) = -1/R3;
G(3,6) = 1;
G(4,4) = -1;
G(4,6) = alpha;
G(5,4) = 1/R4;
G(5,5) = -(R4+R0)/(R4*R0);
G(6,2) = -1;
G(6,3) = 1;

F = zeros(6,1);

% DC Sweep
Vin = linspace(-10,10,1000);
VDCo = zeros(1000,1);
for n = 1 : length(Vin)
    F(1) = Vin(n);
    V = G\F;
    VDCo(n) = V(5);
end

% AC Sweep
F(1) = 1 + 0j; %fix input voltage at 1V AC
VACo = zeros(100,1);
omega = linspace(0,1000,1000);
for n = 1 : length(omega)
    B = G + (1j * omega(n) * C);
    V = B\F;
    VACo(n) = V(5);
end

A = 20*log10(abs(VACo));
semilogx(omega, A);

% Perturbation of C1
N = 10000;
C1p = normrnd(C1,0.05,N,1);
VACp = zeros(N,1);

for n = 1 : length(C1p)
    C(2,1) = C1p(n);
    C(2,2) = -C1p(n);
    B = G + (1j * pi * C);
    V = B\F;
    VACp(n) = V(5);
end

Ap = 20*log10(abs(VACp));

%plot everything
fig = tiledlayout(2,2);
ax_DC = nexttile;
plot(ax_DC,Vin,VDCo);
grid on
xlabel('V_{in} (V)');
ylabel('V_{out} (V)');
title('DC Sweep');

ax_AC = nexttile;
semilogx(ax_AC,omega,A);
grid on
xlabel('\omega (rad/s)');
ylabel('Gain (dB)');
title('AC Sweep');

ax_Cp = nexttile;
histogram(ax_Cp, C1p);
xlabel('C (F)');
title('Perturbed Capacitance');

ax_Gp = nexttile;
histogram(ax_Gp, Ap);
xlabel('Gain (dB)');
title('Perturbed Gain')
