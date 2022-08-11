clc;clear;close all;
id = 192291;

% linearity.m => A <0.9
A = 0.88;
A_step = 0.3;
Ts = 2;
% pocty vzorku nezavisle na Ts
n_ident = 3500; %pro identifikaci z PRBS
n_val = 1500; %pro validaci
n_step = 50; %pro step
n_quiet = 100; %pro oddeleni PRBS a step
upperT = (n_ident+n_val+n_step+n_quiet)*Ts;
upperl = (n_ident+n_val+n_step+n_quiet);
PRBSband = 0.09; %PRBS band
%casovy vektor
t = 0:Ts:upperT;
% vstupni vektor step, quiet zone, PRBS
u = zeros(1,length(t));
u(1:n_step+1) = ones(1,n_step+1)*A_step;
u((n_step+n_quiet+1):(upperl)) = ...
    idinput((n_val+n_ident), 'prbs', [0 PRBSband], [-A A])';
% ziskani vystupu
y = odezva_2021(id, u, t);
y = y';
% rozdeleni zpet dle urceni
%step
t_step = 0:Ts:n_step*Ts;
y_step = y(1:n_step+1);
u_step = u(1:n_step+1);
%ident
y_ident = y((n_step+n_quiet+1):(upperl+1-n_val));
u_ident = u((n_step+n_quiet+1):(upperl+1-n_val));
t_ident = 0:Ts:n_ident*Ts;
%validace
y_val = y((upperl+1-n_val):end);
u_val = u((upperl+1-n_val):end);
t_val = 0:Ts:n_val*Ts;
%
y = y_ident';
u = u_ident';
t = t_ident;

%% obyc MNC
k = 3:length(u);
PHI = [-y(k-1), -y(k-2), u(k-1), u(k-2)];
Y = y(k);
% chytre skrty
for i=length(Y):-1:1
   if Y(i) == 0 || PHI(i,1) == 0
       Y(i) = [];
       PHI(i,:) = [];
   end
end
TH1 = PHI\Y;
y2 = zeros(size(u));
for k = 3:length(u)
      y2(k) = -y2(k-1)*TH1(1)-y2(k-2)*TH1(2)+u(k-1)*TH1(3)+u(k-2)*TH1(4);
end

%% MNC se zpozdenim pozorovanim
d = 4;
k = 4 + d:length(u);
PHI = [-y(k-1), -y(k-2), u(k-1), u(k-2)];
DZ = [-y(k-1-d),-y(k-2-d), u(k-1), u(k-2)];
Y = y(k);
% chytre skrty
for i=length(Y):-1:1
   if Y(i) == 0 || PHI(i,1) == 0 || DZ(i,1) == 0
       Y(i) = [];
       PHI(i,:) = [];
       DZ(i,:) = [];
   end
end
TH2 = (DZ'*PHI)\(DZ'*Y);
y4 = zeros(size(Y));
for k = 3:length(u)
      y4(k) = -y4(k-1)*TH2(1)-y4(k-2)*TH2(2)+u(k-1)*TH2(3)+u(k-2)*TH2(4);
end

%% dodatecny model
y_ivm = zeros(1, length(y));
for i = 3:length(y_ivm)
    y_ivm(i) = -y_ivm(i-1)*TH1(1)-y_ivm(i-2)*TH1(2)+u(i)*TH1(3)+u(i-1)*TH1(4);
end
k = 3:length(u);
Y = y(k);
PHI = [-y(k-1), -y(k-2) u(k-1), u(k-2)];
DZ = [-y_ivm(k-1)', -y_ivm(k-2)', u(k-1), u(k-2)];
% chytre skrty
for i=length(Y):-1:1
   if Y(i) == 0 || PHI(i,1) == 0 || DZ(i,1) == 0
       Y(i) = [];
       PHI(i,:) = [];
       DZ(i,:) = [];
   end
end
TH3 = (DZ'*PHI)\(DZ' * Y);
y5 = zeros(size(u));
for k = 3:length(u)
    y5(k) = -y5(k-1)*TH3(1)-y5(k-2)*TH3(2)+u(k-1)*TH3(3)+u(k-2)*TH3(4);
end

%% ident toolbox
load('192291.mat')
tbx = tf([ident_model.B],[ident_model.A],Ts);
y_val_tbx = lsim(tbx, u_val, t_val);
y_step_tbx = lsim(tbx, u_step, t_step);

%% validace
% MNC
MNC = tf([TH1(3) TH1(4)],[1 TH1(1) TH1(2)],Ts);
y_MNC_val = lsim(MNC, u_val, t_val);
y_MNC_step = lsim(MNC, u_step, t_step);

% MNCSZP
MNCSZP = tf([TH2(3) TH2(4)],[1 TH2(1) TH2(2)],Ts);
y_MNCSZP_val = lsim(MNCSZP, u_val, t_val);
y_MNCSZP_step = lsim(MNCSZP, u_step, t_step);

% MNCSDM
MNCSDM = tf([TH3(3) TH3(4)],[1 TH3(1) TH3(2)],Ts);
y_MNCSDM_val = lsim(MNCSDM, u_val, t_val);
y_MNCSDM_step = lsim(MNCSDM, u_step, t_step);

%% PLOT
% step
figure(1)
plot(t_step, y_step, t_step, y_MNC_step, t_step, y_MNCSZP_step,...
    t_step, y_MNCSDM_step, t_step, y_step_tbx)
title('Odezva na skokovou změnu')
legend('Odezva','MNC','MNCSZP','MNCSDM','Toolbox','Location','southeast')
% identifikace
figure(2)
subplot(4,1,1)
plot(t,u,'k',t,y,'r')
title('Odezva na PRBS')
legend('PRBS', 'Odezva','Location','southeastoutside')

subplot(4,1,2)
plot(t,u,'k',t,y2,'b',t,y,'r')
title('Metoda nejmenších čtverců')
legend('PRBS', 'MNČ','Reálný','Location','southeastoutside')

subplot(4,1,3)
plot(t,u,'k',t,y4,'b',t,y,'r')
title('MNČ se zpožděným pozorováním')
legend('PRBS', 'MNČSZP','Reálný','Location','southeastoutside')

subplot(4,1,4)
plot(t,u,'k',t,y5,'b',t,y,'r')
title('MNČ s dodatečným modelem')
legend('PRBS', 'MNČSDM','Reálný','Location','southeastoutside')
% validace
figure(3)
subplot(5,1,1)
plot(t_val,u_val,'k',t_val,y_val,'r')
title('Odezva na PRBS')
legend('PRBS', 'Odezva','Location','southeastoutside')

subplot(5,1,2)
plot(t_val,u_val,'k',t_val,y_MNC_val,'b',t_val,y_val,'r')
title('Validace: Metoda nejmenších čtverců')
legend('PRBS', 'MNČ','Reálný','Location','southeastoutside')

subplot(5,1,3)
plot(t_val,u_val,'k',t_val,y_MNCSZP_val,'b',t_val,y_val,'r')
title('Validace: MNČ se zpožděným pozorováním')
legend('PRBS', 'MNČSZP','Reálný','Location','southeastoutside')

subplot(5,1,4)
plot(t_val,u_val,'k',t_val,y_MNCSDM_val,'b',t_val,y_val,'r')
title('Validace: MNČ s dodatečným modelem')
legend('PRBS', 'MNČSDM','Reálný','Location','southeastoutside')

subplot(5,1,5)
plot(t_val,u_val,'k',t_val,y_val_tbx,'b',t_val,y_val,'r')
title('Validace: Ident Toolbox')
legend('PRBS', 'Ident TBX','Reálný','Location','southeastoutside')
