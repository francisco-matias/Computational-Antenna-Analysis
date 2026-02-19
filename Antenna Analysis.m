clear;
clf;
close all;

% Constants
u0 = 4*pi*10^-7;          % Vacuum permeability [H/m]
e0 = 8.8542e-12;          % Vacuum permittivity [F/m]

% Given data
V = 1;                    % Input voltage [V]
lambda = 1;               % Wavelength [m]
a = 0.001*lambda;         % Dipole radius [m]
d = 10*a;
k0 = (2*pi)/lambda;       % Wave number [rad/m]
n0 = sqrt(u0/e0);         % Vacuum impedance

% Data for ex1 and ex2
L = 0.47*lambda;          % Dipole length [m]
N = [21, 51, 71];         % Number of steps

Zin = zeros(1,3);         % Input impedance

for j = 1:length(N)

    Ns = N(j);
    delta = L/(2*Ns);     
    zm = zeros(1,Ns);     

    % Matrices
    Zs = zeros(Ns,Ns);    
    Zc = zeros(Ns,Ns);    
    wc = zeros(Ns,1);
    ws = zeros(Ns,1);

    column_zeros = zeros(Ns,1);
    row_zeros = zeros(1,Ns);
    row_u = zeros(1,Ns);

    % Calculate segment midpoints
    for n = 1:Ns
        zm(n) = (n - 0.5)*delta;
    end

    for n = 1:Ns

        % Kernel functions
        Kc = @(z) ((1/(4*pi)) .* ...
            (exp(-1i*k0*sqrt(a^2 + (zm(n)-z).^2)) ./ sqrt(a^2 + (zm(n)-z).^2) + ...
             exp(-1i*k0*sqrt(a^2 + (zm(n)+z).^2)) ./ sqrt(a^2 + (zm(n)+z).^2)));

        Kself = @(z) ((1/(4*pi)) .* ...
            (exp(-1i*k0*sqrt(d^2 + (zm(n)-z).^2)) ./ sqrt(d^2 + (zm(n)-z).^2) + ...
             exp(-1i*k0*sqrt(d^2 + (zm(n)+z).^2)) ./ sqrt(d^2 + (zm(n)+z).^2)));

        for m = 1:Ns

            Zs(n,m) = integral(Kself, zm(m)-delta/2, zm(m)+delta/2, 'ArrayValued', true);
            Zc(n,m) = integral(Kc,   zm(m)-delta/2, zm(m)+delta/2, 'ArrayValued', true);

            wc(m) = -cos(k0*zm(m));
            ws(m) =  sin(k0*abs(zm(m)));

            if m == Ns
                row_u(m) = 1;
            end
        end
    end

    % Constructing matrix 1
    matrix_1 = [Zs Zc wc column_zeros;
                Zc Zs column_zeros wc;
                row_u row_u 0 0;
                row_zeros row_zeros sin(k0*L/2) -sin(k0*L/2)];

    % Constructing matrix 2
    pre_matrix = [ws;
                  column_zeros;
                  0;
                  cos(k0*L/2)];

    matrix_2 = (-1i*V/(2*n0)) * pre_matrix;

    % Solving system of equations
    final_results = linsolve(matrix_1, matrix_2);

    % Extracting constants and currents
    constant1 = final_results(Ns-1,1);
    constant2 = final_results(Ns,1);
    current1 = final_results(1:Ns);
    current2 = final_results(Ns+1:2*Ns);

    % Calculating input impedance
    Zin(j) = V / final_results(1,1);

    % Generating values for plot
    values = delta*(1:Ns);

    % Plotting currents
    figure(1);
    title('Absolute values of current I_1 for different Ns');
    set(gca,'FontSize',18)
    stairs(values, abs(current1)*1e3, 'LineWidth',1.5);
    hold on;

    figure(2);
    title('Absolute values of current I_2 for different Ns');
    set(gca,'FontSize',18)
    stairs(values, abs(current2)*1e3, 'LineWidth',1.5);
    hold on;

end

% Labels and legends
figure(1);
legend('Ns=21','Ns=51','Ns=71');
xlabel('2z/L');
ylabel('|I_1(z)| [mA]');

figure(2);
legend('Ns=21','Ns=51','Ns=71');
xlabel('2z/L');
ylabel('|I_2(z)| [mA]');

disp('Input impedance values:');
disp(Zin);


%constants:
u0 = 4*pi*10^-7; %vacuum permeability [H/m]
e0 = 8.8542*10^-12; %vacuum permitivity [F/m]
%dados do enunciado:
V = 1; %input voltage [V]
lambda = 1;
a = 0.001*lambda;
d = 10*a;
k0 = (2*pi)/lambda;
n0 = sqrt(u0/e0);

% Dados do ex1 e ex2
L = 0.47*lambda;
N = [21, 51, 71];

Zin = zeros(1,3);   % question 2

for j = 1:length(N)   % passar por todos os Ns

    Ns = N(j);
    delta = L/(2*Ns);     % formula (7)
    zm = zeros(1,Ns);     % formula (7)

    % Matrizes necessarias
    Zs = zeros(Ns,Ns);    % formula (12a)
    Zc = zeros(Ns,Ns);    % formula (12b)
    wc = zeros(Ns,1);
    ws = zeros(Ns,1);

    column_zeros = zeros(Ns,1);
    row_zeros = zeros(1,Ns);
    row_u = zeros(1,Ns);

    % Formula (7)
    for n = 1:Ns
        zm(n) = (n-0.5)*delta;
    end

    for n = 1:Ns   % Kc e Kself

        Kc = @(z) (1/(4*pi)) .* ...
            ( exp(-1i*k0*sqrt(a^2 + (zm(n)-z).^2)) ./ sqrt(a^2 + (zm(n)-z).^2) + ...
              exp(-1i*k0*sqrt(a^2 + (zm(n)+z).^2)) ./ sqrt(a^2 + (zm(n)+z).^2) );

        Kself = @(z) (1/(4*pi)) .* ...
            ( exp(-1i*k0*sqrt(d^2 + (zm(n)-z).^2)) ./ sqrt(d^2 + (zm(n)-z).^2) + ...
              exp(-1i*k0*sqrt(d^2 + (zm(n)+z).^2)) ./ sqrt(d^2 + (zm(n)+z).^2) );

        for m = 1:Ns   % preencher matrizes

            Zs(n,m) = integral(Kself, zm(m)-delta/2, zm(m)+delta/2, 'ArrayValued', true);
            Zc(n,m) = integral(Kc,   zm(m)-delta/2, zm(m)+delta/2, 'ArrayValued', true);

            wc(m) = -cos(k0*zm(m));
            ws(m) =  sin(k0*abs(zm(m)));

            if m == Ns
                row_u(m) = 1;
            end
        end
    end

    % Concatenar matrizes → matriz 1 (formula 14)
    matrix_1 = [Zs Zc wc column_zeros;
                Zc Zs column_zeros wc;
                row_u row_u 0 0;
                row_zeros row_zeros sin(k0*L/2) -sin(k0*L/2)];

    % Matriz 2 (formula 14)
    pre_matrix = [ws;
                  column_zeros;
                  0;
                  cos(k0*L/2)];

    matrix_2 = (-1i*V/(2*n0)) * pre_matrix;

    % Sistema de equacoes
    final_results = linsolve(matrix_1, matrix_2);

    % Guardar constantes
    constant1 = final_results(Ns-1,1);
    constant2 = final_results(Ns,1);

    % Correntes
    current1 = final_results(1:Ns);
    current2 = final_results(Ns+1:2*Ns);

    % Impedancia de entrada (question 2)
    Zin(j) = V / final_results(1,1);

    % Valores para grafico
    values = zeros(1,Ns);
    for n = 1:Ns
        values(n) = delta*n;
    end

    % Grafico
    figure(1);
    title('Absolute values of the currents 1 and 2 of a folded dipole antenna');
    set(gca,'FontSize',18)
    stairs(values, abs(current1)*1e3, 'LineWidth',1.5);
    hold on;
    stairs(values, abs(current2)*1e3, 'LineWidth',1.5);
    hold on;

end

legend('Ns=21: |I1|','Ns=21: |I2|', ...
       'Ns=51: |I1|','Ns=51: |I2|', ...
       'Ns=71: |I1|','Ns=71: |I2|');

xlabel('2z/L');
ylabel('|I(z)| [mA]');

disp('Input impedance values (Q2):');
disp(Zin);

%% =========================
%  CONSTANTES
%% =========================
u0 = 4*pi*1e-7;          
e0 = 8.8542e-12;         
n0 = sqrt(u0/e0);        

%% =========================
%  DADOS DO PROBLEMA
%% =========================
V = 1;                    
lambda = 1;               
k0 = 2*pi/lambda;         
a = 0.001*lambda;         
d = 10*a;                 

%% =========================
%  PARAMETROS NUMERICOS
%% =========================
N = 51;                   
N_pontos = 100;           

L_lambda = linspace(0.1,1.1,N_pontos);
L = L_lambda*lambda;

directivity = zeros(1,N_pontos);

%% =========================
%  LOOP EM L
%% =========================
for j = 1:N_pontos
    
    Ns = N;
    delta = L(j)/(2*Ns);
    zm = ((1:Ns) - 0.5)*delta;

    Zs = zeros(Ns,Ns);
    Zc = zeros(Ns,Ns);
    wc = zeros(Ns,1);
    ws = zeros(Ns,1);
    column_zeros = zeros(Ns,1);
    row_zeros = zeros(1,Ns);
    row_u = zeros(1,Ns);

    %% =========================
    %  CALCULO MATRIZES
    %% =========================
    for n = 1:Ns
        
        Kc = @(z) (1/(4*pi)) * ...
            ( exp(-1i*k0*sqrt(a^2 + (zm(n)-z).^2)) ./ sqrt(a^2 + (zm(n)-z).^2) ...
            + exp(-1i*k0*sqrt(a^2 + (zm(n)+z).^2)) ./ sqrt(a^2 + (zm(n)+z).^2) );

        Kself = @(z) (1/(4*pi)) * ...
            ( exp(-1i*k0*sqrt(d^2 + (zm(n)-z).^2)) ./ sqrt(d^2 + (zm(n)-z).^2) ...
            + exp(-1i*k0*sqrt(d^2 + (zm(n)+z).^2)) ./ sqrt(d^2 + (zm(n)+z).^2) );

        for m = 1:Ns
            
            Zs(n,m) = integral(Kself, zm(m)-delta/2, zm(m)+delta/2);
            Zc(n,m) = integral(Kc, zm(m)-delta/2, zm(m)+delta/2);

            wc(m) = cos(k0*zm(m));
            ws(m) = sin(k0*abs(zm(m)));

            if m == Ns
                row_u(m) = 1;
            end
        end
    end

    %% =========================
    %  SISTEMA LINEAR
    %% =========================
    matrix_1 = [ Zs Zc wc column_zeros;
                 Zc Zs column_zeros wc;
                 row_u row_u 0 0;
                 row_zeros row_zeros sin(k0*L(j)/2) -sin(k0*L(j)/2) ];

    rhs = [ws;
           zeros(Ns,1);
           0;
           cos(k0*L(j)/2)];

    matrix_2 = -(1i*V/(2*n0))*rhs;

    solution = matrix_1 \ matrix_2;

    current1 = solution(1:Ns);
    current2 = solution(Ns+1:2*Ns);

    %% =========================
    %  DIRETIVIDADE
    %% =========================
    theta = linspace(0,pi,200);
    he = zeros(size(theta));

    values = (1:Ns)*delta;

    pp_current1 = pchip(values,current1);
    pp_current2 = pchip(values,current2);

    for t = 1:length(theta)

        theta_rad = theta(t);

        f = @(z) ( (ppval(pp_current1,z) + ...
                    ppval(pp_current2,z)) / current1(1) ) ...
                    .* cos(k0*z*cos(theta_rad));

        he(t) = 2*sin(theta_rad) * ...
            integral(f,0,L(j)/2,'RelTol',1e-6,'AbsTol',1e-9);

        he(t) = abs(he(t));
    end

    numerator = 2*max(he).^2;
    denominator = trapz(theta, abs(he).^2 .* sin(theta));

    directivity(j) = numerator / denominator;
end

%% =========================
%  GRAFICO
%% =========================
figure;
plot(L_lambda, directivity, 'r','LineWidth',1.5);
xlabel('L/\lambda');
ylabel('Directivity');
title('Directivity vs Dipole Length');
grid on;
xlim([0.1 1.1]);
set(gca,'FontSize',14);