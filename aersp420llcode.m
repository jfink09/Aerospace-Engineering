% Jason Fink
% AERSP 420 Lifting Line October 28, 2020

%% Clear command window

clear                           ;
clc                             ;
close all                       ;

%% Define knowns

N = 100                         ;   % # of spanwise segments
a = 2*pi                        ;   % 2D lift-curve slope
AR = 10                         ;   % Aspect ratio = b^2/S
clmax = 1.6                     ;   % Max 2D cL when stall
c = 1/AR                        ;   % Chord for rectangular wing
alpha_array_length = 20         ;   % Alpha array length 

%% Define the size of the arrays

alpha = linspace(0,0.396,20)    ;   % Generate list of alphas

ee = zeros(100,1)               ;   % Initializing ee vector
Q = zeros(100,100)              ;   % Initializing Q matrix
W = zeros(100,100)              ;   % Initializing W matrix
CC = zeros(100,100)             ;   % Initializing CC matrix

gamma = zeros(100,1)            ;   % Initializing gamma vector

Cl = zeros(100,1)               ;   % Initializing Cl matrix
cL = zeros(20,1)                ;   % Initializing cL vector
%% Computing the matrices

for i = 1:1:100                      % Nest for loops to assign elements i = [1,N] j = [1,N]
    
    ee(i) = 1                    ;   % Vector of ones

    
    for j = 1:1:100
        
        Q(i,j) = 1/(1-4*(i-j)^2) ;   % Finding Q
        
        if i == j
            
            I(i,j) = 1           ;   % Identity
            CC(i,j) = c          ;   % Finding C
  
        end
    end 
end
%% More definitions
                        
c_vec = c*ee                     ;   % Vector of chord values
AR_vec = AR*ee                   ;   % Vector of 1/c
CC_diag = diag(c_vec)            ;   % Diagonal chord matrix
AR_diag = diag(AR_vec)           ;   % Diagonal 1/c matrix
%% Finding W

W = ((a*N)/(2*pi))*Q             ;    % Finding W matrix

%% Find CLalpha 

alpha_min = 0                               ;    % Minimum alpha to check
alpha_max = 20*pi/180                       ;    % Maximum alpha to check in rdians
k = 40                                      ;    % Number of steps in alpha sweep
alpha_list = linspace(alpha_min,alpha_max,k);    % Generage list of alphas
cL = zeros(k,1)                             ;    % Initialize list of 3D cLs

        
stallFrac = 0                                     ;  % Initialize stallFrac = 0                           
cLmaxFound = 0                                    ;  % Initialize tracker of if cLmax is found

for i = 1:1:k                                        % Loop to iterate angles of attack
    stallCount = 0                                ;  % Initialize a count of the number of spanwise stations that stalled
    alpha = alpha_list(i)                         ;  % Pick alpha from alpha_list
    gamma = linsolve((AR_diag+W),(0.5*a*alpha*ee));  % Calculate gamma using linsolve function to avoid taking the inverse
    cl = 2*AR_diag*gamma                          ;  % Caclulate local lift coefficient

    for j = 1:1:N                                    % Sweep across spanto check locations that stalled
        
        if cl(j) > clmax                                 
            cl(j) = 0                             ;  % Set cl to zero if greater than clmax
            stallCount = stallCount + 1           ;  % Count locations that stalled
            stallFrac = stallCount/N              ;  % Calculate fraction of span that stalled
            
        end
        
    end
    gamma = 0.5*c*cl                          ;  % Calculate new vorticity distribution accounting for stalled portions of the wing
    cL(i) = (2/(ee'*c*ee))*ee'*gamma          ;  % Calculate 3D cL for the entire wing
        
    if (stallFrac > 0) && (cLmaxFound == 0)      % Check if wing stalled
        cLmaxFound = 1                        ;  % Note that cLmax was found
        cLmax = cL(i)                         ;  % Note what that cL is
    end
end

cLalpha = ((ee'*linsolve((AR_diag+W),ee))/(ee'*c*ee))*a;  % Calculate 3D lift curve slope and avoid taking the inverse 

% Print required calculations 
fprintf('Calculations: \n')
fprintf('cLalpha = %.3f\n', cLalpha)                                    
fprintf('cLmax (using predicted cL) = %.3f\n',max(cL))
fprintf('cLmax (any wing section stall) = %.3f\n',cLmax)    
%% plot CL vs alpha

figure(1)

ax = axes                                                         ;   % define axes                                   
ax.ColorOrder = [0 0 0]                                           ;   % Make plot black

hold on

title('Lift coefficient cL with respect to Angle of Attack (deg)');   % Title
xlabel('Angle of Attack (deg)')                                   ;   % Label x-axis               
ylabel('Lift Coefficient cL')                                     ;   % Label y-axis

plot(alpha_list*(180/pi),cL(:,1),'linewidth',3)                       % Plot lift coefficient vs. angle of attack in degrees
grid on                                                               % Major grid
grid minor                                                            % Minor grid

hold off


%% cLalpha with respect to N plot
% Angle of attack = 1/10

  d = zeros(901,1)                                                        ;   % Initialize the number of spanwise segments matrix as 'd'
  K = 1                                                                   ;   % Set a constant equal to 1
 
  for w = 100:1:1000                                                          % Nest for loops to assign elements w = [100,1000] i = [1,n]
     
      N = w                                                               ;   % Set the spanwise segments matrix equal to w matrix              
    
      for i = 1:1:N                                                           % i = [1,100]
    
          ee(i) = 1                                                       ;   % Vector of ones
     
        for j = 1:1:N                                                         % j = [1,100]
        
            Q(i,j) = 1/(1-4*(i-j)^2)                                      ;   % Define Q matrix
        
            if i == j                                                         % If i equals j set the C matrix equal to the chord
            
                I(i,j) = 1                                                ;   % Define the identity matrix
                CC(i,j) = c                                               ;   % Defne C matrix
            
            end
      
        end
       
    end
    
    W = ((a*N)/(2*pi))*Q                                                  ;   % Define the W matrix
    cLalpha(K) = (1/trace(CC))*ee'*linsolve((linsolve(CC,I) + W),I)*ee*a  ;   % Define the cLalpha matrix 
    d(K) = w                                                              ;   % Define n equal to w for each k
    
    K = K + 1                                                             ;   % Counter

  end
  
%% plot cLalpha vs N

figure(2)

ax = axes                                                                            ;  % define axes
ax.ColorOrder = [0 0 0]                                                              ;  % Make plot black

hold on

title('Lift curve slope with respect to spanwise discretization N (alpha = 1/10)')   ;  % Title
xlabel('N')                                                                          ;  % Label x-axis
ylabel('cLalpha')                                                                    ;  % Label y-axis

plot(d,cLalpha,'linewidth',3)                                                        ;  % Plot N with lift coefficient

grid on                                                                                 % Major grid
grid minor                                                                              % Minor grid

hold off
