% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%
% Nom du fichier d'output a analyser (modifiez selon vos besoins)
filename = 'output.out'; 

% Chargement des donnees
data = load(filename);
%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = ''; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2_2020.exe'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base MODIFIER SELON VOS BESOINS

nsimul = 20; % Nombre de simulations a faire

%nsteps = round(logspace(2,4,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS
nsteps = round(logspace(1,4,nsimul))

paramstr = 'nsteps'; % Nom du parametre a scanner  MODIFIER SELON VOS BESOINS
param = nsteps; % Valeurs du parametre a scanner  MODIFIER SELON VOS BESOINS

%% Simulations %% 
%%%%%%%%%%%%%%%%%
% Lance une serie de simulations (= executions du code C++)
% Normalement, on ne devrait pas avoir besoin de modifier cette partie

output = cell(1, nsimul); % Tableau de cellules contenant les noms des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i)), '.out']
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i})
    disp(cmd)
    system(cmd);
  
end

%% Analyse %%
%%%%%%%%%%%%%
% Ici, on aimerait faire une etude de convergence: erreur fonction de dt, sur diagramme log-log.
% A MODIFIER ET COMPLETER SELON VOS BESOINS

%----------------------------Exercie 2.3.a.i------------------------------
%NE PAS MODIFIER 

%Erreur max sur la position

error_pos_max = zeros(1, nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    xend = data(end,2); % Extraire le x final
    zend = data(end,3); % Extraire le z final
    x_th = ((1.6726e-27/(1.6022e-19*5))*4e5*cos(1.6022e-19*5*1.e-8/1.6726e-27)-(1.6726e-27*4e5/(1.6022e-19*5))); % TODO: Entrer la vraie solution analytique a tfin
    z_th = ((4e5*1.6726e-27/(1.6022e-19*5))*sin(1.6022e-19*5*1.e-8/1.6726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    
    error_pos_max(i) = max(abs(xend-x_th), abs(zend-z_th)); % Maximum de l'erreur sur la position
end


lw=1; fs=16;
figure('Name', [filename ': Convergence numérique de l''erreur sur la position'])
plot(dt, error_pos_max, 'b+-','linewidth',lw)
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Delta t')
ylabel('Convergence de l''erreur de la position')
grid on
 
 x = x_th %Pour écrire dans la fenêtre (sans ;)
 z = z_th


%Erreur max sur la vitesse
error_vit_max = zeros(1, nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    v_xend = data(end,4); % Extraire le x final
    v_zend = data(end,5); % Extraire le z final
    v_x_th = (-4e5*sin(1.6022e-19*5*1.e-8/1.62726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    v_z_th = (4e5*cos(1.6022e-19*5*1.e-8/1.62726e-27)); % TODO: Entrer la vraie solution analytique a tfin

    error_vit_max = max(abs(v_xend-v_x_th), abs(v_zend-v_z_th)); % Maximum de l'erreur sur la position
end


figure('Name', [filename ': Convergence numérique de l''erreur sur la vitesse'])
plot(dt, error_vit_max, 'b+-', 'linewidth',lw)
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Delta t')
ylabel('Convergence de l''erreur sur la vitesse')
grid on

v_x = v_x_th
v_z = v_z_th

%Erreur sur l'Energie mécanique 

delta_E_mec = zeros(1, nsimul);
w_c_dt = zeros(1,nsimul) 
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    E_mec_t = data(end,6); % Extraire E_mec(t)
    E_mec_0 = (1.0/2.0)*1.62726e-27*4e5*4e5;
    delta_E_mec(i) = E_mec_t - E_mec_0; % erreur sur la position finale
    w_c_dt(i) = (1.6022e-19 *5 )/1.62726e-27 * dt(i)%Pour l'intégrateur BB
end

figure('Name', [filename ': Convergence numérique de l''erreur sur la vitesse'])
plot(dt, delta_E_mec, 'b+-', 'linewidth',lw)
%hold on  -> Pour Boris Buneman
%plot(dt, w_c_dt, 'r+-', 'linewidth', lw)
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Delta t')
ylabel('\Delta E_{mec}')
grid on


%-----------------Fin de l'exo


