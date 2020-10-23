% Nom du fichier d'output a analyser (modifiez selon vos besoins)
filename = 'output.out'; 

% Chargement des donnees
data = load(filename);

% Extraction des quantites d'interet - Figures de bases 
% (Le code c++ ecrit t, x(t), z(t), v_x(t), v_z(t), E_mec(t), M_magn(t)  ligne par ligne, 
%  une ligne par pas de temps)
t = data(:,1); 
x = data(:,2);
z = data(:,3)
v_x = data(:,4);
v_z = data(:,5)
Emec = data(:,6);
M_magn = data(:,7)

% nombre de pas de temps effectués:
nsteps = length(t)
% longueur du pas de temps:
dt = t(2)-t(1)

% Figures
% line width and font size (utile pour la lisibilité des figures dans le
% rapport)
lw=2; fs=16; 
figure('Name', [filename ': x(t)'])
plot(t, x, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('x [m]')
grid on

lw=2; fs=16; 
figure('Name', [filename ': z(t)'])
plot(t, z, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('Name', [filename ': v_x(t)'])
plot(t, v_x, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('v_x [m/s]')
grid on

figure('Name', [filename ': v_z(t)'])
plot(t, v_z, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('v_z [m/s]')
grid on

lw=0.25;
figure('Name', [filename ': (x,z)'])
plot(x, z, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('x [m]')
ylabel('z [m]')
grid on

lw=0.25;
figure('Name', [filename ': (v_x,v_z)'])
plot(v_x, v_z, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('v_x [m/s]')
ylabel('v_z [m/s]')
grid on

lw=0.25;
figure('Name', [filename ': (x,v_x)'])
plot(x, v_x, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('x [m]')
ylabel('v_x [m/s]')
grid on

lw=0.25;
figure('Name', [filename ': (z,v_z)'])
plot(z, v_z, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('z [m]')
ylabel('v_z [m/s]')
grid on

figure('Name', [filename ': Emec(t)'])
plot(t, Emec, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('Emec [J]')
grid on

figure('Name', [filename ': M_magn(t)'])
plot(t, M_magn, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('M_magn [A*m^2]')
grid on


%Graph de Delta_E_mec

dt = t(2)-t(1);
E_mec_0 = data(1:1, 6);
Delta_E_mec = Emec - E_mec_0;
figure('Name', [filename ': D_E_mec(t)'])

plot(t, Delta_E_mec, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('\Delta E_{mec} [J]')
grid on

