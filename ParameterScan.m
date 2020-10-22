% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

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
error = zeros(1,nsimul);
error_x = zeros(1, nsimul);
error_z = zeros(1, nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    xend = data(end,2); % Extraire le x final
    zend = data(end,3); % Extraire le z final
    x_th = ((1.6726e-27/(1.6022e-19*5))*4e5*cos(1.6022e-19*5*1.e-8/1.6726e-27)-(1.6726e-27*4e5/(1.6022e-19*5))); % TODO: Entrer la vraie solution analytique a tfin
    z_th = ((4e5*1.6726e-27/(1.6022e-19*5))*sin(1.6022e-19*5*1.e-8/1.6726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    error(i) = sqrt((xend-x_th)^2+(zend-z_th)^2); % erreur sur la position finale
    error_x(i) = abs(xend-x_th); % Calcul de l'erreur locale en x
    error_z(i) = abs(zend-z_th); % Calcul de l'erreur locale en z
    
    erreur_pos_max = max(error_x(i), error_z(i)) % Maximum de l'erreur sur la position
end




figure
loglog(dt, error_x, 'k+',"Color", "blue")
hold on
loglog(dt, error_z, 'k+', "Color", "magenta")
xlabel('\Delta t')
ylabel('Erreur locale sur la position (x(b), z(m))')
grid on
 
 x = x_th %Pour écrire dans la fenêtre (sans ;)
 z = z_th


%Erreur max sur la vitesse
error = zeros(1,nsimul);
error_v_x = zeros(1, nsimul);
error_v_z = zeros(1, nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    v_xend = data(end,4); % Extraire le x final
    v_zend = data(end,5); % Extraire le z final
    v_x_th = (-4e5*sin(1.6022e-19*5*1.e-8/1.62726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    v_z_th = (4e5*cos(1.6022e-19*5*1.e-8/1.62726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    error(i) = sqrt((v_xend-v_x_th)^2+(v_zend-v_z_th)^2); % erreur sur la position finale
    error_v_x(i) = abs(v_xend-v_x_th); %Calcul de l'erreur locale en x
    error_v_z(i) = abs(v_zend-v_z_th); %Calcul de l'erreur locale en z

    erreur_vit_max = max(error_v_x(i), error_v_z(i)) % Maximum de l'erreur sur la position
end


figure
loglog(dt, error_v_x, 'k+',"Color", "blue")
hold on
loglog(dt, error_v_z, 'k+', "Color", "magenta")
xlabel('\Delta t')
ylabel('Erreur locale de la vitesse (v_x(b), v_z(m))')
grid on

v_x = v_x_th
v_z = v_z_th


%-----------------Fin de l'exo


