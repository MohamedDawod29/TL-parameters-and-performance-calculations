clc;
clear All;

%TASK 1
% Entering the conductor parameters
fprintf('HINT:ALL DISTANCES,DIAMETERS AND LENGTHS MUST BE POSITIVE!\n');
rho = input('Enter conductor resistivity: ');
L = input('Enter conductor length in km: ');
L = L * 1000;
d = input('Enter conductor diameter in cm: ');
d = d/100;
% Calculating conductor radius and cross-sectional area
r = d/2;
A = pi*(r)^2;

% Choosing symmetrical or unsymmetrical spacing
spacing = menu('Choose if the TL is symmetrical spacing or unsymmetrical spacing','Symmetrical','Unsymmetrical');

if spacing == 1
    % The system is symmetrical
    D = input('Enter distance between phases: ');
    while D <= 0
       fprintf('The distance must be greater than 0\n');
       D = input('Enter distance between phases: '); 
    end   
    
    % Calculating resistance, inductance, and capacitance per phase for symmetrical
    R_DC = rho*L/A;
    R_AC = 1.1*R_DC;
    L1 = 2e-7*log(D/(0.7788*r));
    C_ph = 2*pi*8.854e-12/log(D/r);
else
    % The system is unsymmetrical
    Dab = input('Enter distance between phases a and b: ');
     while Dab <= 0
       fprintf('The distance must be greater than 0\n');
       Dab = input('Enter distance between phases a and b: '); 
     end 
    
    Dbc = input('Enter distance between phases b and c: ');
     while Dbc <= 0
       fprintf('The distance must be greater than 0\n');
       Dbc = input('Enter distance between phases b and c: '); 
    end   
    
    Dca = input('Enter distance between phases c and a: ');
     while Dca <= 0
       fprintf('The distance must be greater than 0\n');
       Dca = input('Enter distance between phases c and a: ');
     end   
    
    Deq = power(Dab*Dbc*Dca, 1/3);
    
    % Calculating resistance, inductance, and capacitance per phase for unsymmetrical
    R_DC = rho*L/A;
    R_AC = 1.1*R_DC;
     
     L1 = (2e-7*log(Deq/(0.7788*r)));
     C_ph = (2*pi*8.854e-12/log(Deq/r));
end

fprintf('the value of the DC resistance = %d ohm\n',R_DC); 
fprintf('the value of the AC resistance = %d ohm\n', R_AC); 
fprintf('the value of the inductance per phase = %d H \n',L1); 
fprintf('the value of the capacitance per phase = %d F\n',C_ph);

%TASK 2
f=50;
XL=2*pi*f*L1*L;            
Y=2*pi*f*C_ph*L*1i;  
Z=R_AC+1i*XL;

if(L<80000)
    A=1;
    B=Z;
    C=0;
    D=1;
elseif(L>=80000 && L<=250000)
model = menu('Choose the model you want?(pi or T)','pi-model','T-model');
    switch model
        case 1
            A=1+(Z.*Y)./2;
            B=Z;
            C=Y.*(1+(Z.*Y)./4);
            D=1+(Z.*Y)./2;
        case 2
            A=1+(Z.*Y)./2;
            B=Z.*(1+(Z.*Y)./4);
            C=Y;
            D=1+(Z.*Y)./2;             
     end   
else
    fprintf("ADCD parameters can't be calculated\n");
end

A
B
C
D

%TASK 3
VrA = input('please input recieved voltage amplitude (in KVs) \nNote it will be considered phase voltage\n');
VrA = VrA * 1000;
Vr = ((VrA) / sqrt(3)) * exp(0);

fprintf('please select which case you would like to carry out\n'); 
caseselect = input ('input 1 for case 1 (varying power) and input 2 for case 2 (varying power factor) \n');
while true
    if (caseselect ~= 1) && (caseselect ~= 2)
        caseselect = input ('error. please pick 1 or 2 \n');
        continue
    end
    break
end

if caseselect == 1
    %case 1
    Pr = 0: 1000: 100000;
    Sr = Pr / (0.8) * exp(1i * acos(0.8));
    Ir = conj(Sr / (3 * Vr));
    
    Vs = A * Vr + B * Ir;
    Is = C * Vr + D * Ir;
    
    Vreg = (abs(Vs / A) - abs (Vr)) / abs (Vr);
    figure
    plot(Pr, Vreg*100);
    xlabel('Active power Received');
    ylabel('Voltage Regulation');
    title('Voltage Regulations vs Varying Active Power Recieved') 
    Ps= 3 * abs(Vs) .* abs(Is) .* cos(angle (Vs .* conj(Is)));
    eff=Pr ./ Ps;
    figure
    plot (Pr, eff*100);
    xlabel('Active power');
    ylabel('Efficiency');
    title('Efficiency vs Varying Active Power Recieved') 
else
    %case 2 which is varying the pf instead of the precieved

    Pr=100000;

    %Lagging Power Factors
    pf = linspace (0.3,1,100); 
    Sr1 = Pr ./ (pf) .* exp(1i * acos(pf));
    Ir1 = conj(Sr1 / (3 * Vr));

    Vs1 = A * Vr + B * Ir1;
    Is1 = C * Vr + D * Ir1;
    
    Ps1 = 3 * abs(Vs1) .* abs(Is1) .* cos(angle(Vs1 .* conj(Is1)));

    Vreg1 = ((abs(Vs1 / A) - abs(Vr)) / abs(Vr)) * 100;
    eff1 = (Pr ./ Ps1) * 100;
    
    %Leading Power Factors
    Sr2 = Pr ./ (pf) .* exp(1i * -acos(pf));    
    Ir2 = conj(Sr2 / (3 * Vr));
    
    Vs2 = A * Vr + B * Ir2;
    Is2 = C * Vr + D * Ir2;

    Vreg2 = ((abs(Vs2 / A) - abs(Vr)) / abs(Vr)) * 100;
    Ps2 = 3 * abs(Vs2) .* abs (Is2) .* cos (angle (Vs2 .* conj(Is2)));
    eff2 = (Pr ./ Ps2) * 100;
    
    subplot(2,2,1)
    plot (pf, Vreg1);
    xlabel('Power Factor Lagging');
    ylabel('Voltage Regulation (%)');
    title('Voltage Regulation vs Varying Power Factor (Lag to UPF)');
    
    subplot(2,2,2);
    plot(pf, eff1);
    xlabel('Power Factor Lagging');
    ylabel('Efficiency (%)');
    title('efficiency with varying power factors (Lag to UPF)');

    subplot(2,2,3);
    plot(pf, Vreg2);
    xlabel('Power Factor Leading');
    ylabel('Voltage Regulation (%)');
    title('Voltage Regulation vs Varying Power Factors (Lead to UPF)');
    
    subplot(2,2,4);
    plot(pf, eff2);
    xlabel('Power Factor Leading');
    ylabel('Efficiency (%)');
    title('Efficiency vs Varying Power Factor (Lead to UPF)')

end