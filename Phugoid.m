%% Clear workspace

clc                     ;
clear                   ;
close all               ;

%% Define variables and constants

Iy = 1249               ;
e = 0.6                 ;
epsilon = 0.413         ;
at = 3.66               ;
chord = 5.25            ;
aw = 5.05               ;
b = 30                  ;
Sw = b*chord            ;
weight = 2400           ;
m = weight/32           ;
g = 32                  ;
f = 5                   ;
lt = 13.96-1            ;
St = 10*2.5             ;
eta = 1                 ;
S = St+Sw               ;
AR = Sw/b               ;

%% Calculations

Cl_alpha = aw+eta*at*(St/S)*(1-epsilon)                                         ;

hnw = 0.25                                                                      ;
hl = lt/chord                                                                   ;

hn = (hnw + hl*eta*(St/S)*(at/aw)*(1-epsilon))/(1+hl*(St/S)*(at/aw)*(1-epsilon));

distToCG = 86/12                                                                ;
distToLE = 78.4/12                                                              ;

h = (distToCG - distToLE)/chord                                                 ;

StabMargin = (hn-h)                                                             ;

Cm_alpha = (h-hn)*Cl_alpha                                                      ;

fprintf('Cl_alpha = %f\n',Cl_alpha)                                             ;
fprintf('hn = %f\n',hn)                                                         ;
fprintf('h = %f\n',h)                                                           ;
fprintf('StabMargin = %f\n',StabMargin)                                         ;
fprintf('Cm_alpha = %f\n',Cm_alpha)                                             ;


%% Stability derivatives

U0 = 147                                                ;
theta0 = 0                                              ;
rho = 0.00237                                           ;

drag = 0.5*rho*U0.^2*f                                  ;
fprintf('drag = %f\n',drag)                             ;


T0 = drag                                               ;
fprintf('T0 = %f\n',T0)                                 ;

Xu = (-1/m)*((T0/U0)+rho*U0*f)                          ;
fprintf('Xu = %f\n',Xu)                                 ;

Cl0 = weight/(0.5*rho*U0.^2*S)                          ;
fprintf('Cl0 = %f\n',Cl0)                               ;

pi = 3.1415926535                                       ;

X_alpha = (-1/m)*((rho*U0.^2*Sw*Cl0*Cl_alpha)/(pi*AR*e));
fprintf('X_alpha = %f\n',X_alpha)                       ;

Zu = (-1/m)*rho*U0*Sw*Cl0                               ;
fprintf('Zu = %f\n',Zu)                                 ;

Z_alpha = (-1/m)*0.5*rho*U0.^2*Sw*Cl_alpha              ;
fprintf('Z_alpha = %f\n',Z_alpha)                       ;

Zq = -eta*(1/(2*m))*rho*U0*at*St*lt                     ;
fprintf('Zq = %f\n',Zq)                                 ;

M_alpha = (1/Iy)*0.5*rho*U0.^2*S*chord*Cm_alpha         ;
fprintf('M_alpha = %f\n',M_alpha)                       ;

Mq = (m/Iy)*lt*Zq                                       ;
fprintf('Mq = %f\n',Mq)                                 ;

Z_alphadot = -eta*(1/(2*m))*rho*U0*at*St*lt*epsilon     ;
fprintf('Z_alphadot = %f\n',Z_alphadot)                 ;

M_alphadot = (m/Iy)*lt*Z_alphadot                       ;
fprintf('M_alphadot = %f\n\n',M_alphadot)               ;

%% Matrix form

MatrixForm = ["U0", U0; 
              "theta0", theta0; 
              "Iy", Iy; 
              "e", e; 
              "epsilon", epsilon; 
              "at", at; 
              "chord", chord; 
              "Cl_alpha", Cl_alpha; 
              "Sw", Sw; 
              "weight", weight; 
              "m", m; 
              "g", g; 
              "f", f; 
              "rho", rho; 
              "lt", lt; 
              "St", St; 
              "eta", eta; 
              "f", f; 
              "M_alpha", M_alpha];
fprintf('Variables in matrix form: \n');
disp(MatrixForm)


%% State Space

AA = [Xu, X_alpha, -g, 0; (Zu/(U0-Z_alphadot)), (Z_alpha/(U0-Z_alphadot)), 0, ((Zq+U0)/(U0-Z_alphadot)); 
      0, 0, 0, 1; 
      ((M_alphadot*Zu)/(U0-Z_alphadot)), M_alpha+((M_alphadot*Z_alpha)/(U0-Z_alphadot)), 0, Mq+M_alphadot*((Zq+U0)/(U0-Z_alphadot))];

fprintf('AA state space matrix form: \n');
disp(AA)

%% Eigenvalues

eigenvalues = eig(AA)                                               ;
fprintf('Eigenvalues: \n')                                          ;
disp(eigenvalues)

% Create variable for each eigenvalue to be added to a table
e1 = eigenvalues(1,:)                                               ;
e2 = eigenvalues(2,:)                                               ;
e3 = eigenvalues(3,:)                                               ;
e4 = eigenvalues(4,:)                                               ;

% Show eigenvalues in table form 
Table = table({'e1';'e2';'e3';'e4'},...
    [e1;e2;e3;e4],...
    'VariableNames',{'Eigenvalue Labels' 'Eigenvalues'},...
    'RowNames',{'e1' 'e2' 'e3' 'e4'})                                     

%% Find the real part of the eigenvalues

real_e1 = real(e1);
real_e2 = real(e2);
real_e3 = real(e3);
real_e4 = real(e4);

%% Load csv data and make plots

T = readtable('phugoidData.xlsx','Range','A1:E942')                   ;

time = table2array(T(:,1))                                            ;
A5 = table2array(T(:,5))                                              ;

figure(1)
fplot(@(t) real(exp(e1*t)), [0 3],'Color','#EDB120','lineWidth',2)
hold on 
grid on
grid minor
fplot(@(t) real(exp(e2*t)), [0 3],'Color','#EDB120', 'lineWidth',2)
title('Amplitude vs. Time (s)')
xlabel('Time (s)')                                                    ;
ylabel('Amplitude')                                                   ;
hold off

figure(2)
fplot(@(t) real(exp(e3*t)), [0 200],'Color','#0072BD','lineWidth',2)
hold on
grid on
grid minor
title('Amplitude vs. Time (s)')
xlabel('Time (s)')                                                    ;
ylabel('Amplitude')                                                   ;
hold off

figure(3)
fplot(@(t) real(exp(e1*t)), [0 20],'Color','#0072BD','lineWidth',2)
hold on
grid on
grid minor
fplot(@(t) real(exp(e3*t)), [0 20],'Color','#EDB120','lineWidth',2)
title('Amplitude vs. Time (s)')
xlabel('Time (s)')                                                    ;
ylabel('Amplitude')                                                   ;

%% Define McCormick values
g2 = 32                ; 
theta02 = 0            ; 
U02 = 147              ; 
X_alpha2 = -27.64      ; 
Z_alpha2 = -258        ; 
Zq2 = -2.655           ; 
M_alphadot2 = -0.7149  ; 
Xu2 = -0.07151         ; 
Zu2 = -0.439           ; 
M_alpha2 = -17.34      ; 
Mq2 = -1.927           ;
Z_alphadot2 = -1.1375  ;

%% Comparison with McCormick with table
NewTable = ["Xu", Xu, Xu2; 
              "Zu", Zu, Zu2; 
              "X_alpha", X_alpha, X_alpha2; 
              "Z_alpha", Z_alpha, Z_alpha2; 
              "M_alpha", M_alpha, M_alpha2; 
              "Zq", Zq, Zq2; 
              "Mq", Mq, Mq2; 
              "Z_alphadot", Z_alphadot, Z_alphadot2; 
              "M_alphadot", M_alphadot, M_alphadot2];

fprintf('Comparison with McCormick: \n');
disp(NewTable)

%% Plots

figure(4)
hold on
xlabel('Time (s)')                                                                              ;
ylabel("g's")                                                                                   ;
title("g's vs. Time (s)")
grid on
grid minor 
plot(time,A5,'.')
hold off

figure(5)
fplot(@(t) -1-0.6*(exp(real(e3)*t)*cos(imag(e3)*t)), [0 90],'Color','#0072BD','lineWidth',2)
hold on
grid on
grid minor
title("g's vs. Time (s)")
xlabel('Time (s)')                                                                              ;
ylabel("g's")                                                                                   ;

figure(6)
fplot(@(t) -1-0.6*(exp(real(e3)*t)*cos(imag(e3)*t)), [0 90],'Color','#0072BD','lineWidth',2)
hold on
grid on
grid minor
plot(time,A5,'.','Color','#0072BD')
title("g's vs. Time (s)")
xlabel('Time (s)')                                                                              ;
ylabel("g's")                                                                                   ;

