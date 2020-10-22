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
ylabel('Emec [s]')
grid on

figure('Name', [filename ': M_magn(t)'])
plot(t, M_magn, '-','linewidth',lw)
set(gca,'fontsize',fs)
xlabel('t [s]')
ylabel('M_magn [A*m^2]')
grid on


%-------------------------Exercice2.3----------------------------------

%-------------------------Partie a-------------------------------------

%----------------------------i)-----------------------------------------

%Convergence de l'erreur de la position -------- Euler explicite

nsimul  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
error_pos = [0.0014, 
             9.7108e-04,
             5.8968e-04,
             3.8472e-04,
             2.5421e-04,
             1.6935e-04,
             1.1468e-04,
             7.8801e-05,
             5.3915e-05,
             3.7003e-05,
             2.5599e-05,
             1.7684e-05,
             1.2259e-05,
             8.5042e-06,
             5.9026e-06,
             4.0989e-06,
             2.8475e-06,
             1.9786e-06,
             1.3750e-06,
             9.5563e-07
    ];

lw=2; fs=16;
figure('Name', [filename ': Erreur sur la position selon Euler explicite'])
plot(1./nsimul, error_pos, 'k+-','linewidth',lw, "Color", "blue")
set(gca,'fontsize',fs)
xlabel('1/N simulation')
ylabel('error position max')
grid on

%Convergence de l'erreur de la vitesse----------Euler explicite
nsimul  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20];
error_pos = [6.9847e+05,
            4.7275e+05,
            2.9008e+05,
            1.9191e+05,
            1.2940e+05,
            8.8760e+04,
            6.2577e+04,
            5.0867e+04,
            5.1260e+04,
            5.1639e+04,
            5.1947e+04,
            5.2185e+04,
            5.2360e+04,
            5.2486e+04,
            5.2576e+04,
            5.2640e+04,
            5.2685e+04,
            5.2717e+04,
            5.2739e+04,
            5.2754e+04];

lw=2; fs=16;
figure('Name', [filename ': Erreur sur la vitesse selon Euler explicite'])
plot(1./nsimul, error_pos, 'k+-','linewidth',lw, "Color", "blue")
set(gca,'fontsize',fs)
xlabel('1/N simulation')
ylabel('error vitesse max')
grid on

%-----------------------------ii)---------------------------------------
% Non-conservation de l'énergie mécanique 
E_mec_i = (1.0/2.0 * 1.62726e-27*4e5*4e5);
E_mec_t = [7.75697638819292e-007,
            7.76455503400251e-007,
            7.76990057899759e-007,
            7.7736157050707e-007,
            7.77614331728026e-007,
            7.77780741836893e-007 ]
delta_E_mec = E_mec_t-E_mec_0;
plot(nsteps_num, abs(delta_E_mec),'k+-', "Color", "blue")
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('N_{steps}')
ylabel('delta E_mec(t)')
grid on

