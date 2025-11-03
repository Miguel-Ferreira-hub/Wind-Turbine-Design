%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA EXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
filename='CLD25.txt';
filename_reference='500kW Reference Turbine Blade profile Test.txt';
data=readtable('CLD25.txt');
data_reference=readtable('500kW Reference Turbine Blade profile Test.txt');
angle=data.Var1;
CL=data.Var2;
CD=data.Var3;
r=data_reference.radialLengthR_m_;
c=data_reference.chordLengthC_m_;
beta=data_reference.twistAngleBeta_degree_;
areadata=readtable('Cross-sectional area.txt');
spars=readtable('Spars.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEGMENT SPLIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sf=6.0976; %18MW
%sf=4.2927; %9MW
r_new=r.*sf;
c_new=c.*sf;


combined_r=linspace(r_new(1), r_new(17),19)';
combined_c = interp1(r.*sf, c.*sf, combined_r, 'spline');
combined_beta = interp1(r.*sf, beta, combined_r, 'spline');

new_data = table(combined_r, combined_beta, combined_c, ...
    'VariableNames', {'RadialLength_m', 'TwistAngle_deg', 'ChordLength_m'});

bladearea=areadata.Cross_sectionalArea;
bredth=spars.SparCapBredth_B_m_(1:19);
depth=spars.WebDepth_2d_m_(1:19);
thickness=spars.ChosenThickness_m_(1:19);
spararea=(2.*bredth.*thickness) + (2.*(depth-2.*thickness).*thickness);
rhoCFRP=1580;
rhoGFRP=1850;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALL FUNCTIONS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vo=11;
omega=0.6160;
rpm=omega*60/(2*pi);
rotor_frequency=rpm/60;
bladepassing_frequency=3*rpm/60;
fprintf('RPM: %e\n', rpm);
fprintf('Rotor Frequency (1P - Hz): %e\n', rotor_frequency);
fprintf('Blade Passing Frequency (3P - Hz): %e\n', bladepassing_frequency);

%BEM_f(combined_r, combined_c, combined_beta, 3, 11, 0, true);

%BEM_f(r, c, beta, 3, 9.5, 0, true); %Test 0.5MW

%Performance(5.5, 30, 20, combined_r, combined_c, combined_beta, 3, 0);

%Pitchangle_Performance(5.5, 25, 100, combined_r, combined_c, combined_beta, 3);

%[CoM, TotalMoment, totalmass]=A(combined_r, bladearea, spararea, rhoGFRP, rhoCFRP);

[f]=D(combined_r, bladearea, spararea, rhoGFRP, rhoCFRP, 3, 6, 80e-3, 7850, 165, 200e9);

%disp(rotor_frequency/f);

%envelope();

tabletower=table(6, 80, 7850, 165, 200, 11800, 20600*3, 'VariableNames', {'Average Tower Diameter (m)', 'Tower Thickness (mm)', 'Density (kg/m^3)', 'Tower Height (m)', 'Youngs Modulus (GPa)', 'Tower Mass (kg)', 'Hub Mass (kg)'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [angle, CL, CD]=coefficients(i)
 %s=[33,37,17, 41, 46, 29, 44, 15, 21, 8, 19, 20, 31, 36, 48, 32, 28,25, 12
 
 %s=[33,37,17, 41, 46, 29, 44, 21, 15, 8, 19, 20, 31, 36, 48, 32, 28, 12];
 s=[33,37,17, 41, 46, 29, 44, 21, 15, 8, 19, 20, 31, 36, 48, 32, 2, 28, 12]; %actual one
 %s=[33, 17, 16, 46, 7, 24, 44, 6, 21, 8, 19, 18, 31, 36, 48, 32, 28, 25]; %2
 %s=[7, 24, 44, 35, 6, 21, 18, 10, 34, 25, 9, 26, 5, 1, 23, 7, 27 48]; %1
 % el bueno da 14 
  % s=[33,37,17, 41, 46, 29, 44, 21, 15, 8, 19, 20, 31, 36, 48, 32, 28,12, 47 ];
%s=[33,37,44,6,21,34,8,26,23,18,20,31,40,36,48,32,2,49,28,25,12];

    filename = sprintf('CLD%02d.txt', s(i));  % Select file based on segment index i
    
    % Read data, skipping the first row (header)
    data = readmatrix(filename, 'NumHeaderLines', 1);
    
    % Extract data
    angle = data(:, 1);
    CL = data(:, 2);
    CD = data(:, 3);
end

function [aa, ar]=induction_factors(v0, omega, r, theta, sigma, angle, CD, CL, n)
    aa=zeros(n:1);
    ar=zeros(n:1);
    xi=(omega*r)/v0;
    aa(1)=0;
    ar(1)=0;
    for i=2:n
        phi=atan(((1-aa(i-1))/(1+ar(i-1)))/xi);
        AOA=phi*180/pi-theta;
        CL_int=interp1(angle,CL,AOA,'spline');
        CD_int=interp1(angle,CD,AOA,'spline');
        CN=CL_int*cos(phi)+CD_int*sin(phi);
        CR=CL_int*sin(phi)-CD_int*cos(phi);
        aa(i)=1/(4*sin(phi)*sin(phi)/(sigma*CN)+1);
        ar(i)=1/(4*sin(phi)*cos(phi)/(sigma*CR)-1);
    end
    fprintf('aa: %e\n', aa(end));
    fprintf('ar: %e\n', ar(end));
end

function [aa, ar, F, K, CN, CR, phi]=induction_factors_with_tiplossfactor(v0, omega, r, theta, B, c, angle, CD, CL, R, n, Display)
    aa=zeros(n,1);
    ar=zeros(n,1);
    F=zeros(n,1);
    K=zeros(n,1);
    xi=(omega*r)/v0;
    aa(1)=0;
    ar(1)=0;
    F(1)=0;
    K(1)=0;
    ac=1/3;
    s=c*B/(2*pi*r);
    CN=zeros(n,1);
    CR=zeros(n,1);
    CN(1)=0;
    CR(1)=0;
    CL_int=zeros(n,1);
    CD_int=zeros(n,1);
    CL_int(1)=0;
    CD_int(1)=0;
    phi=zeros(n,1);
    phi(1)=0;
    for i=2:n
        phi(i)=atan(((1-aa(i-1))/(1+ar(i-1)))/xi);
        %fprintf('i=%d, phi=%f deg\n', i, phi(i)*180/pi);
        AOA=(phi(i)*180/pi)-theta;
        CL_int(i)=interp1(angle,CL,AOA,'spline');
        CD_int(i)=interp1(angle,CD,AOA,'spline');
        CN(i)=CL_int(i)*cos(phi(i))+CD_int(i)*sin(phi(i));
        CR(i)=CL_int(i)*sin(phi(i))-CD_int(i)*cos(phi(i));
        F(i) = 2/pi*acos(exp(-B/2*(R-r)./r./sin(phi(i))));
        K(i) = 4*F(i).*(sin(phi(i))).^2./(s.*CN(i));
        aa(i) = 1./(K(i)+1);
        if aa(i) > ac
           aa(i) = (1-(sqrt(1+4/K(i)*((1-ac)/(1-2*ac))^2)-1)*(K(i)*(1-2*ac)/2));
        end
        ar(i)=1/((4*F(i)*sin(phi(i))*cos(phi(i))/(s*CR(i)))-1);

        if ~isreal(phi(i)) || ~isreal(CL_int(i)) || ~isreal(CD_int(i)) || ~isreal(CN(i)) || ~isreal(CR(i)) || ...
   ~isreal(F(i)) || ~isreal(K(i)) || ~isreal(aa(i)) || ~isreal(ar(i))

        phi(i)=NaN;
        CL_int(i)=NaN;
        CD_int(i)=NaN;
        CN(i)=NaN;
        CR(i)=NaN;
        F(i) = NaN;
        K(i) = NaN;
        aa(i) = NaN;
        ar(i)=NaN;
        end
        
    end

    if Display==true
        fprintf('aa: %e\n', aa(end));
        fprintf('ar: %e\n', ar(end));
        fprintf('F: %e\n', F(end));
        fprintf('K: %e\n', K(end));
        fprintf('CN: %e\n', CN(end));
        fprintf('CR: %e\n', CR(end));
        fprintf('Phi: %e\n', phi(end));
    end

end

function [fr, fn, T, tau, Power, Cp, Cplimit, shell, axialinductionfactor, angularinductionfactor]=BEM_f(r, c, beta, B, v0, pitch, Display)
n=10;
omega=0.6160;  %18MW
%omega=0.8750;  %9MW
P=101000;
R=287;
Temp=287.5;
rho=1.225;
Mflap=zeros(1,length(r)-1);
Medge=zeros(1,length(r)-1);
shellthickness=zeros(1,length(r));
axialinductionfactor=zeros(1,length(r));
angularinductionfactor=zeros(1,length(r));
tiploss=zeros(1,length(r));
phiarray=zeros(1,length(r));
alphaarray=zeros(1,length(r));
fn=zeros(1,length(r));
fr=zeros(1,length(r));
T=zeros(1,length(r)-1);
tau=zeros(1,length(r)-1);
ri=zeros(1, length(r)-1);
xiprime=zeros(1, length(r)-1);
Rend=20.5*6.0976;
%Rend=20.5*4.2927;
    for i=1:length(r)

        [angle, CL, CD] = coefficients(i);  

        theta=beta(i)+pitch;
        [aa, ar, F, K, CN, CR, phi] = induction_factors_with_tiplossfactor(v0, omega, r(i), theta, B, c(i), angle, CD, CL, Rend, n, false);
        axialinductionfactor(i)=aa(end);
        angularinductionfactor(i)=ar(end);
        tiploss(i)=F(end);
        phiarray(i)=phi(end);
        alphaarray(i)=(phi(end).*180./pi)-theta;
        vrel=(1-aa(end))./sin(phi(end))*v0;
        fn(i)=0.5*CN(end)*rho*vrel^2*c(i);
        fr(i)=0.5*CR(end)*rho*vrel^2*c(i);
        %fprintf('i=%d, CR=%f, fr=%f\n', i, CR(end), fr(i));

        %Print Statements
        if Display==true
        fprintf('Element: %0.f\n', i);
        fprintf('aa: %e\n', aa(end));
        fprintf('ar: %e\n', ar(end));
        fprintf('F: %e\n', F(end));
        fprintf('K: %e\n', K(end));
        fprintf('CN: %e\n', CN(end));
        fprintf('CR: %e\n', CR(end));
        fprintf('Relative Velocity: %e\n', vrel);
        fprintf('fr: %e\n', fr(i));
        fprintf('fn: %e\n', fn(i));
        end
    end
    for i=1:length(r)-1
    %Average Length Along Turbine
        ri(i)=((r(i)+r(i+1))/2);
        xiprime(i)=ri(i)*omega/v0;

    %Forces
        %Normal Force
        T(i)=(0.5)*((fn(i+1)+fn(i))*(r(i+1)-r(i)));
        %Rotational Force
        tau(i)=(1/6)*(((fr(i+1)+fr(i))*(r(i+1)^2-r(i)^2))+(fr(i+1)*r(i+1)^2-fr(i)*r(i)^2)-((fr(i+1)-fr(i))*(r(i+1)*r(i))));
    
    %Moments
        %Flap-wise Moments
        Mflap(i)=(0.5)*fn(i)*((r(i+1)^2)-(r(i)^2));
        %Edge-wise Moments
        Medge(i)=(0.5)*fr(i)*((r(i+1)^2)-(r(i)^2));
    end

     Mflap_total=flip(cumsum(flip(Mflap)));
     Medge_total=flip(cumsum(flip(Medge)));
     
     for i=1:length(r)-1
     rinitial=50e-3;
     shellthickness(i)=rinitial*(Mflap_total(i)/max(Mflap_total));
     if shellthickness(i)<10e-3
         shellthickness(i)=10e-3;
     end

     end
     shellthickness(19)=10e-3;
     segment=linspace(1, 19, 19);
     shell=table(segment(:), shellthickness(:), 'VariableNames', {'Segment', 'Shell Thickness'});
     
     %disp('flap');
     %disp(Mflap_total(:));
     %disp('edge');
     %disp(Medge_total(:));
%Power and Total Forces
TotalT=sum(T);
Totaltau=sum(tau);
Power=B*omega*(Totaltau);
Cp=Power/(0.5*rho*(v0^3)*pi*(Rend^2));
Cplimit=Cp*(27/16);

if Display==true
fprintf('Total Normal Force [kN]: %e\n', TotalT/1000);
fprintf('Total Rotational Force [kN]: %e\n', Totaltau/1000);
fprintf('Power [W]: %e\n', Power);
fprintf('Power Coefficient (Cp): %e\n', Cp);
fprintf('Power Coefficient Without Correction (Cp): %e\n', Cplimit);
fprintf('Total Length: %e\n', Rend);

%fn and fr plot
figure('Name', 'Forces');
plot(r, fn./1000, 'o-', 'LineWidth', 1);
hold on
plot(r, fr./1000, 'o-', 'LineWidth', 1);
text(80, 9.5, '$f_{N,i}$', 'Interpreter', 'latex', 'Color', [0 0.447 0.741], 'FontSize', 14);
text(80, 2.3, '$f_{R,i}$', 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098], 'FontSize', 14);
xlabel('$r_i$, m', 'Interpreter', 'latex');
ylabel('$f_N,i, f_R,i$, kN/m', 'Interpreter', 'latex');
title('Normal and Rotational Force Coefficients', 'Interpreter', 'latex');
hold off

%plot of Total Forces
figure('Name', 'Total Forces');
plot(xiprime, T./1000, 'o-', 'LineWidth', 1);
hold on
plot(xiprime, tau./1000, 'o-', 'LineWidth', 1);
text(4.5, 100, '$T_i$', 'Interpreter', 'latex', 'Color', [0 0.447 0.741], 'FontSize', 14);
text(4.5, 720, '$\tau_i$', 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098], 'FontSize', 14);
xlabel('$\bar{\xi} = \Omega \bar{r_i} / V_0$', 'Interpreter', 'latex');
ylabel('$\tau_i$, kNm; $T_i$, kN', 'Interpreter', 'latex');
title('Normal and Rotational Forces', 'Interpreter', 'latex');
hold off

%plot of Chord Length and Beta distribution
figure('Name', 'Chord Length and Beta Distribution');
plot(r, c, 'o-', 'LineWidth', 1);
hold on
plot(r, beta, 'o-', 'LineWidth', 1);
text(80, 8, '$c$', 'Interpreter', 'latex', 'Color', [0 0.447 0.741], 'FontSize', 14);
text(80, 4.2, '$\beta$', 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098], 'FontSize', 14);
xlabel('$r_i$, m', 'Interpreter', 'latex');
ylabel('$c_i$, m, $\beta_i, ^ \circ$', 'Interpreter', 'latex');
title('Chord Length and Twist Angle Distribution', 'Interpreter', 'latex');
hold off

%plot of Flap-wise Bending Moments
figure('Name', 'Flapwise Bending Moments');
plot(ri, Mflap_total./1000, 'o-', 'LineWidth', 1);
hold on
plot(ri, Mflap./1000, 'o-', 'LineWidth', 1);
xlim([29 122]);
text(90, 55100, '$M_y(\bar{r_i})$', 'Interpreter', 'latex', 'Color', [0 0.447 0.741], 'FontSize', 14);
text(75, 15000, '$M_y,i,$', 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098], 'FontSize', 14);
xlabel('$\bar{r_i}$, m', 'Interpreter', 'latex');
ylabel('$M_y,i, M_y(\bar{r_i})$, kNm', 'Interpreter', 'latex');
title('Flapwise Bending Moments', 'Interpreter', 'latex');
hold off

%plot of Edge-wise Bending Moments
figure('Name', 'Flapwise Bending Moments');
plot(ri, Medge_total./1000, 'o-', 'LineWidth', 1);
hold on
plot(ri, Medge./1000, 'o-', 'LineWidth', 1);
text(84, 20000, '$M_z(\bar{r_i})$', 'Interpreter', 'latex', 'Color', [0 0.447 0.741], 'FontSize', 14);
text(70, 5000, '$M_z,i,$', 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098], 'FontSize', 14);
xlim([29, 122]);
xlabel('$\bar{r_i}$, m', 'Interpreter', 'latex');
ylabel('$M_z,i, M_z(\bar{r_i})$, kNm', 'Interpreter', 'latex');
title('Edgewise Bending Moments', 'Interpreter', 'latex');
hold off

%plot of Edge-wise Bending Moments
figure('Name', 'Shell Thickness');
plot(r, shellthickness.*1000, 'o-', 'LineWidth', 1);
ylabel('Shell Thickness, mm', 'Interpreter', 'latex');
xlabel('$r_i, m$', 'Interpreter', 'latex');
title('Shell Thickness against Radial Distance', 'Interpreter', 'latex');

%Verification Code
figure('Name', 'Axial Induction Factor');
plot(r, axialinductionfactor, 'o-', 'LineWidth', 1);
ylabel('$a$', 'Interpreter', 'latex');
xlabel('$r_i, m$', 'Interpreter', 'latex');
title('Axial Induction Factors', 'Interpreter', 'latex');

figure('Name', 'Angular Induction Factor');
plot(r, angularinductionfactor, 'o-', 'LineWidth', 1);
ylabel('$a^{\prime}$', 'Interpreter', 'latex');
xlabel('$r_i, m$', 'Interpreter', 'latex');
title('Angular Induction Factors', 'Interpreter', 'latex');

figure('Name', 'Tip Loss Factor');
plot(r, tiploss, 'o-', 'LineWidth', 1);
ylabel('$F$', 'Interpreter', 'latex');
xlabel('$r_i, m$', 'Interpreter', 'latex');
title('Tip Loss Factors', 'Interpreter', 'latex');

figure('Name', 'Phi');
plot(r, phiarray.*180./pi, 'o-', 'LineWidth', 1);
hold on
plot(r, alphaarray, 'o-', 'LineWidth', 1);
hold off
text(64, 12, '$\phi$', 'Interpreter', 'latex', 'Color', [0 0.447 0.741], 'FontSize', 14);
text(60, 6.7, '$\alpha$', 'Interpreter', 'latex', 'Color', [0.85 0.325 0.098], 'FontSize', 14);
ylabel('$\phi,\ \alpha = \phi - \beta,\ ^\circ$', 'Interpreter', 'latex');
xlabel('$r_i, m$', 'Interpreter', 'latex');
title('Angle of Attack (AOA) and Phi', 'Interpreter', 'latex');

end
end

function Performance(vmin, vmax, n, r, c, beta, B, pitch)
         R=20.5*6.0976;
         omega=0.6160;
         %omega=0.84;
         v=linspace(vmin, vmax, n);
         vline=linspace(0, 11, 10);
         Px=linspace(0, 18.5, 10);
         vx=zeros(1, length(Px))+11;
         Pline=zeros(1, length(vline))+18.5;
         Powerarray=zeros(1, n);
         Cparray=zeros(1, n);
         Cplimitarray=zeros(1, n);
         lambda=zeros(1, n);
         axialInductionMatrix = zeros(n, length(r));
         angularInductionMatrix = zeros(n, length(r));
         for i=1:n
             [~, ~, ~, ~, Power, Cp, Cplimit, ~, axialinductionfactor, angularinductionfactor]=BEM_f(r, c, beta, B, v(i), pitch, false);
             Powerarray(i)=Power;
             Cparray(i)=Cp;
             Cplimitarray(i)=Cplimit;
             lambda(i)=(omega*R)/v(i);
             disp(axialinductionfactor);
             %disp(angularinductionfactor);
             if istable(axialinductionfactor)
             a = table2array(axialinductionfactor(:,1)); % fix here
             else
             a = axialinductionfactor;
             end
             if istable(angularinductionfactor)
             ar = table2array(angularinductionfactor(:,1)); % fix here
             else
             ar = angularinductionfactor;
             end

             axialInductionMatrix(i, :) = a;  % transpose to match row shape
             angularInductionMatrix(i, :) = ar;

         end
     
figure;

% Create grid
[ri_grid, V0_grid] = meshgrid(r, v);  

% Axial Induction Factor Plot
mesh(ri_grid, V0_grid, axialInductionMatrix);
set(gca, 'YDir', 'reverse');
xlabel('$r_i$, m', 'Interpreter','latex');
ylabel('$V_0$, m/s', 'Interpreter','latex');
zlabel('$a$', 'Interpreter', 'latex');
title('Axial Induction Factor $a$', 'Interpreter', 'latex');

% Angular Induction Factor Plot
figure;
mesh(ri_grid, V0_grid, angularInductionMatrix);
set(gca, 'XDir', 'reverse');
xlabel('$r_i$, m', 'Interpreter','latex');
ylabel('$V_0$, m/s', 'Interpreter','latex');
zlabel('$a''$', 'Interpreter', 'latex');
title('Angular Induction Factor $a''$', 'Interpreter', 'latex');


figure('Name', 'Performance analysis as nominal wind varies');
plot(v, Powerarray./1000000, 'o-', 'LineWidth', 1);
hold on
plot(vline, Pline, 'r', 'LineWidth', 1);
plot(vx, Px, 'r', 'LineWidth', 1);
%plot(v, dashedline, '--', 'Color', 'k');
xlabel('$V_0$, m/s', 'Interpreter', 'latex');
ylabel('P, MW', 'Interpreter', 'latex');
title('Power at Different Wind Speeds', 'Interpreter', 'latex');
%xlim([5 20]); 
%ylim([0 600]);
hold off

figure('Name', 'Performance analysis as nominal wind varies');
plot(v, Cparray, 'o-', 'LineWidth', 1);
xlabel('$V_0$, m/s', 'Interpreter', 'latex');
ylabel('Cp', 'Interpreter', 'latex');
title('Power Coefficient as Wind Varies', 'Interpreter', 'latex');
%xlim([5 20]); 
%ylim([0 0.5]);

figure('Name', 'Performance analysis as nominal wind varies');
plot(lambda, Cplimitarray, 'o-', 'LineWidth', 1);
xlabel('$\lambda = \Omega R / V_0$', 'Interpreter', 'latex');
ylabel('$Cp \times 27 / 16$', 'Interpreter', 'latex');
title('Power Coefficient as Wind Varies', 'Interpreter', 'latex');
%ylim([0 1]);
end

function Pitchangle_Performance(vmin, vmax, n, r, c, beta, B)
         v=linspace(vmin, vmax, n);
         vline=linspace(0, vmax, n);
         pitch=[0, 5, 8, 12, 16, 20];
         %dashedline=zeros(1, n)+460;
         dashedline=zeros(1, n)+18.5;
         %line=zeros(n, 1)+469.113;
         %liney=linspace(0, 800, n);
         envelope=zeros(1, 34);
         varray=zeros(1, 34);
         Powerarray1=zeros(1, n);
         Powerarray2=zeros(1, n);
         Powerarray3=zeros(1, n);
         Powerarray4=zeros(1, n);
         Powerarray5=zeros(1, n);
         Powerarray6=zeros(1, n);
         Cparray1=zeros(1, n);
         Cparray2=zeros(1, n);
         Cparray3=zeros(1, n);
         Cparray4=zeros(1, n);
         Cparray5=zeros(1, n);
         Cparray6=zeros(1, n);
         Tarray1=zeros(1, n);
         Tarray2=zeros(1, n);
         Tarray3=zeros(1, n);
         Tarray4=zeros(1, n);
         Tarray5=zeros(1, n);
         Tarray6=zeros(1, n);
         for i=1:n
             [~, ~, T1, ~, Power1, Cp1, ~, ~, ~, ~]=BEM_f(r, c, beta, B, v(i), pitch(1), false);
             Powerarray1(i)=Power1;
             Cparray1(i)=Cp1;
             Tarray1(i)=sum(T1);
             if Powerarray1(i) < 0 || Cparray1(i) < 0 || Tarray1(i) < 0
                 Powerarray1(i)=NaN;
                 Cparray1(i)=NaN;
                 Tarray1(i)=NaN;
             end
             [~, ~, T2, ~, Power2, Cp2, ~, ~, ~, ~]=BEM_f(r, c, beta, B, v(i), pitch(2), false);
             Powerarray2(i)=Power2;
             Cparray2(i)=Cp2;
             Tarray2(i)=sum(T2);
             if Powerarray2(i) < 0 || Cparray2(i) < 0 || Tarray2(i) < 0
                 Powerarray2(i)=NaN;
                 Cparray2(i)=NaN;
                 Tarray2(i)=NaN;
             end
             [~, ~, T3, ~, Power3, Cp3, ~, ~, ~, ~]=BEM_f(r, c, beta, B, v(i), pitch(3), false);
             Powerarray3(i)=Power3;
             Cparray3(i)=Cp3;
             Tarray3(i)=sum(T3);
             if Powerarray3(i) < 0 || Cparray3(i) < 0 || Tarray3(i) < 0
                 Powerarray3(i)=NaN;
                 Cparray3(i)=NaN;
                 Tarray3(i)=NaN;
             end
             [~, ~, T4, ~, Power4, Cp4, ~, ~, ~, ~]=BEM_f(r, c, beta, B, v(i), pitch(4), false);
             Powerarray4(i)=Power4;
             Cparray4(i)=Cp4;
             Tarray4(i)=sum(T4);
             %if (v(i)<7)
                %Cparray4(i)=NaN;
             %end
             if Powerarray4(i) < 0 || Cparray4(i) < 0 || Tarray4(i) < 0
                 Powerarray4(i)=NaN;
                 Cparray4(i)=NaN;
                 Tarray4(i)=NaN;
             end
             [~, ~, T5, ~, Power5, Cp5, ~, ~, ~, ~]=BEM_f(r, c, beta, B, v(i), pitch(5), false);
             Powerarray5(i)=Power5;
             Cparray5(i)=Cp5;
             Tarray5(i)=sum(T5);
             if Powerarray5(i) < 0 || Cparray5(i) < 0 || Tarray5(i) < 0
                 Powerarray5(i)=NaN;
                 Cparray5(i)=NaN;
                 Tarray5(i)=NaN;
             end
             [~, ~, T6, ~, Power6, Cp6, ~, ~, ~, ~]=BEM_f(r, c, beta, B, v(i), pitch(6), false);
             Powerarray6(i)=Power6;
             Cparray6(i)=Cp6;
             Tarray6(i)=sum(T6);
             if Powerarray6(i) < 0 || Cparray6(i) < 0 || Tarray6(i) < 0                 
                 Powerarray6(i)=NaN;
                 Cparray6(i)=NaN;
                 Tarray6(i)=NaN;
             end
             
         end

         %Performance Envelope
         envelope(1:29)=Powerarray1(1:29)./1e6;
         varray(1:29)=v(1:29);
         envelope(30:34)=[18.4, 18.5, 18.59, 18.5, 18.5];
         varray(30:34)=[11.83, 13.09, 15.33, 18.07, 21.26];

figure('Name', 'Performance analysis as nominal wind varies');
plot(v, Powerarray1./1000000, 'o-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1);
text(4.5, 85, '$\theta_p = 0^\circ$', 'Interpreter', 'latex', 'Color', [0, 0.4470, 0.7410], 'FontSize', 14);
hold on
plot(vline, dashedline, '--', 'Color', 'k');
plot(v, Powerarray2./1000000, 'o-', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1);
text(4.5, 79, '$\theta_p = 5^\circ$', 'Interpreter', 'latex', 'Color', [0.4940, 0.1840, 0.5560], 'FontSize', 14);
plot(v, Powerarray3./1000000, 'o-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
text(4.5, 73, '$\theta_p = 8^\circ$', 'Interpreter', 'latex', 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 14);
plot(v, Powerarray4./1000000, 'o-', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1);
text(4.5, 67, '$\theta_p = 12^\circ$', 'Interpreter', 'latex', 'Color', [0.4660, 0.6740, 0.1880], 'FontSize', 14);
plot(v, Powerarray5./1000000, 'o-', 'Color', [0, 0.2470, 0.5410], 'LineWidth', 1);
text(4.5, 61, '$\theta_p = 16^\circ$', 'Interpreter', 'latex', 'Color', [0, 0.2470, 0.5410], 'FontSize', 14);
plot(v, Powerarray6./1000000, 'o-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1);
text(4.5, 55, '$\theta_p = 20^\circ$', 'Interpreter', 'latex', 'Color', [0.6350, 0.0780, 0.1840], 'FontSize', 14);
xlabel('$V_0$, m/s', 'Interpreter', 'latex');
ylabel('P, MW', 'Interpreter', 'latex');
title('Power as Wind Speed Varies at Different Pitch', 'Interpreter', 'latex');
xlim([4 22]); 
%ylim([0 870]);
hold off

figure('Name', 'Performance Envelope');
plot(varray, envelope, 'o-', 'LineWidth', 1);
xlabel('$V_0$, m/s', 'Interpreter', 'latex');
ylabel('$Power, MW$', 'Interpreter', 'latex');
text(6.5, 11, '$\theta_p = 0^\circ$', 'Interpreter', 'latex', 'FontSize', 10);
text(11, 17, '$\theta_p = 5^\circ$', 'Interpreter', 'latex', 'FontSize', 10);
text(13, 20, '$\theta_p = 8^\circ$', 'Interpreter', 'latex', 'FontSize', 10);
text(15, 17, '$\theta_p = 12^\circ$', 'Interpreter', 'latex', 'FontSize', 10);
text(18, 20, '$\theta_p = 16^\circ$', 'Interpreter', 'latex', 'FontSize', 10);
text(20, 17, '$\theta_p = 20^\circ$', 'Interpreter', 'latex', 'FontSize', 10);
title('$Control Envelope for Nominal Power Output$', 'Interpreter', 'latex');
ylim([0 22]);

figure('Name', 'Performance analysis as nominal wind varies');
plot(v, Cparray1, 'o-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1);
text(6.2, 0.45, '$\theta_p = 0^\circ$', 'Interpreter', 'latex', 'Color', [0, 0.4470, 0.7410], 'FontSize', 14);
hold on
plot(v, Cparray2, 'o-', 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1);
text(8.5, 0.40, '$\theta_p = 5^\circ$', 'Interpreter', 'latex', 'Color', [0.4940, 0.1840, 0.5560], 'FontSize', 14);
plot(v, Cparray3, 'o-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
text(10, 0.31, '$\theta_p = 8^\circ$', 'Interpreter', 'latex', 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 14);
plot(v, Cparray4, 'o-', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1);
text(11, 0.15, '$\theta_p = 12^\circ$', 'Interpreter', 'latex', 'Color', [0.4660, 0.6740, 0.1880], 'FontSize', 14);
plot(v, Cparray5, 'o-', 'Color', [0, 0.2470, 0.5410], 'LineWidth', 1);
text(14, 0.09, '$\theta_p = 16^\circ$', 'Interpreter', 'latex', 'Color', [0, 0.2470, 0.5410], 'FontSize', 14);
plot(v, Cparray6, 'o-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1);
text(17, 0.05, '$\theta_p = 20^\circ$', 'Interpreter', 'latex', 'Color', [0.6350, 0.0780, 0.1840], 'FontSize', 14);
xlabel('$V_0$, m/s', 'Interpreter', 'latex');
ylabel('Cp', 'Interpreter', 'latex');
title('Power Coefficient as Wind Speed Varies at Different Pitch', 'Interpreter', 'latex');
%xlim([5 20]); 
%ylim([0 0.5]);
hold off

figure('Name', 'Performance analysis as nominal wind varies');
plot(Powerarray1./1000000, Tarray1./1000, 'o-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1);
text(10, 1000, '$\theta_p = 0^\circ$', 'Interpreter', 'latex', 'Color', [0, 0.4470, 0.7410], 'FontSize', 14);
hold on
%plot(line, liney, 'Color', 'k');
plot(Powerarray3./1000000, Tarray3./1000, 'o-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
text(22, 500, '$\theta_p = 8^\circ$', 'Interpreter', 'latex', 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 14);
plot(Powerarray6./1000000, Tarray6./1000, 'o-', 'Color', [0.9911, 0.8039, 0.4902], 'LineWidth', 1);
text(18, 180, '$\theta_p = 20^\circ$', 'Interpreter', 'latex', 'Color', [0.9911, 0.8039, 0.4902], 'FontSize', 14);
%plot(v, Tarray4./1000, 'o-', 'LineWidth', 1);
%plot(v, Tarray5./1000, 'o-', 'LineWidth', 1);
%plot(v, Tarray6./1000, 'o-', 'LineWidth', 1);
xlabel('P, MW', 'Interpreter', 'latex');
ylabel('T, kN', 'Interpreter', 'latex');
title('Normal Force as Wind Speed Varies at Different Pitch', 'Interpreter', 'latex');
%xlim([0 800]); 
%ylim([0 20]);
hold off

fprintf('Maximum power for zero pitch angle (kN): %e\n', max(Powerarray1./1000));
end

function [CoM, TotalMoment, totalmass]=A(r, bladearea, spararea, rhoblade, rhospar)
        R=20.5*6.0976;
        g=9.81;
        omega=0.6160;
        %totalarea=spararea+bladearea;
        m1=zeros(1, length(r)-1);
        m2=zeros(1, length(r)-1);
        mtotal=zeros(1, length(r)-1);
        d=zeros(1, length(r)-1);
        massmoment=zeros(1, length(r)-1);
        ri=zeros(1, length(r)-1);
        Momentblade=zeros(1, length(r)-1);
        Momentspar=zeros(1, length(r)-1);
        TotalMoment=zeros(1, length(r)-1);
        Gravitationalforce=zeros(1, length(r)-1);
        Centrifugalforce=zeros(1, length(r)-1);
        for i=1:length(r)-1
            ri(i)=((r(i)+r(i+1))/2);
            d(i)=r(i+1)-r(i);
            m1(i)=bladearea(i).*d(i).*rhoblade; %blade
            m2(i)=spararea(i).*d(i).*rhospar; %spar
            mtotal(i)=m1(i)+m2(i);
            Gravitationalforce(i)=mtotal(i)*g;
            Centrifugalforce(i)=mtotal(i)*omega^2*ri(i);
            massmoment(i)=ri(i)*mtotal(i);
            Momentblade(i)=0.5*bladearea(i)*g*rhoblade*(R-ri(i)-27.439)^2;
            Momentspar(i)=0.5*spararea(i)*g*rhospar*(R-ri(i)-27.439)^2;
            TotalMoment(i)=Momentblade(i)+Momentspar(i);
        end
        CoM=sum(massmoment)/sum(mtotal);
        %Gravitationalforce=0.5*omega^2*sum(mtotal)*R^2;
        %Gravitationalforce=sum(mtotal)*g;
        %Centrifugalforce=sum(mtotal)*omega^2*(CoM-r(1));
        fprintf('Centre of Mass from root: %e\n', CoM);
        fprintf('Centre of Mass (CoM) from blade start: %e\n', CoM-r(1));
        fprintf('Total Mass (Kg): %e\n', sum(mtotal));
        fprintf('Gravitational Load: %e\n', Gravitationalforce);
        fprintf('Centrifugal Load: %e\n', Centrifugalforce);
        TotalMoment=TotalMoment(:);
        segment=linspace(1, 18, 18);
        info=table(segment(:), d(:), ri(:), TotalMoment(:), m1(:), m2(:), mtotal(:), Centrifugalforce(:), Gravitationalforce(:), 'VariableNames', {'Segment', 'Length (m)', 'Half-Way Point (m)', 'Gravitational Moment (Edgewise - Nm)', 'Blade Mass (Kg)', 'Spar Mass (Kg)', 'Total Mass (Kg)', 'Centrifugal Force', 'Gravitational Force'});
        disp(info);
        totalmass=sum(mtotal);
        figure('Name', 'Gravitational Loading');
        plot(ri, TotalMoment./1000, 'o-', 'DisplayName', 'Gravitational Bending Moment' ,'LineWidth', 1);
        xlabel('$Radial Distance from Root (r_i), m$', 'Interpreter', 'latex');
        ylabel('Bending Moment, kN', 'Interpreter', 'latex');
        title('Gravitational Bending Moments (Edgewise)', 'Interpreter', 'latex');

        figure('Name', 'Forces');
        plot(ri, Centrifugalforce./1e6, 'o-', 'DisplayName', '$Magnitude of Centrifugal Force$');
        hold on
        plot(ri, Gravitationalforce./1e6, 'o-', 'DisplayName', '$Magnitude of Gravitational Force$');
        plot(ri, (Centrifugalforce+Gravitationalforce)./1e6, 'o-', 'DisplayName', '$Magnitude of Resultant Force$');
        hold off
        xlabel('$Radial Distance from Root (r_i), m$', 'Interpreter', 'latex');
        ylabel('$Force, MN$', 'Interpreter', 'latex');
        legend('location', 'north', 'Interpreter', 'latex', 'FontSize', 14);
        ylim([0 0.8]);

        figure('Name', 'Mass Distribution');
        plot(ri, m1, 'o-', 'DisplayName', '$Blade Mass$');
        hold on
        plot(ri, m2, 'o-', 'DisplayName', '$Spar Mass$');
        plot(ri, mtotal, 'o-', 'DisplayName', '$Total Mass$');
        xlabel('$\bar{r_i}$, m', 'Interpreter', 'latex');
        ylabel('$Mass, kg$', 'Interpreter', 'latex');
        legend('location', 'north', 'Interpreter', 'latex', 'FontSize', 14);
end

function f=D(r, bladearea, spararea, rhoGFRP, rhoCFRP, B, D, t, rhotower, L, E)
    [~, ~, totalmass]=A(r, bladearea, spararea, rhoGFRP, rhoCFRP);
    Massturbine=(B*totalmass);
    Masstower=(pi*rhotower*D*t);
    disp(Masstower);
    I=(pi/8)*D^3*t;
    f=(1/(2*pi)*sqrt((3*E*I)/((0.24*Masstower + Massturbine)*L^3)));
    fprintf('Natural Frequency of Tower (Hz): %e\n', f);
end

%function vibration=D(r, bladearea, spararea, rhoGFRP, rhoCFRP)
    %[CoM1, ~, totalmass1]=A(r(1:6), bladearea(1:6), spararea(1:6), rhoGFRP, rhoCFRP);
    %[CoM2, ~, totalmass2]=A(r(7:12), bladearea(7:12), spararea(7:12), rhoGFRP, rhoCFRP);
    %[CoM3, ~, totalmass3]=A(r(13:19), bladearea(13:19), spararea(13:19), rhoGFRP, rhoCFRP);
%end

function envelope()
    x1=linspace(0, 11, 10);
    x2=[11.83, 13.09, 15.33, 18.07, 21.26];
    x=zeros(1, 15);
    x(1:10)=x1;
    x(11:15)=x2;
    y=zeros(1, 15);
    y(1:10)=0;
    y(11:15)=[5, 8, 12, 16, 20];
    figure('Name', 'Envelope');
    plot(x, y, 'o-', 'LineWidth', 1);
    ylabel('$Pitch Angle, \circ$', 'Interpreter', 'latex');
    xlabel('$Wind Speed, m/s$', 'Interpreter', 'latex');
    title('Pitch Angle Envelope', 'Interpreter', 'latex');
    ylim([0 25]);
end


