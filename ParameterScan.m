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

nsimul = 4; % Nombre de simulations a faire

nsteps = round(logspace(2,4,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS

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
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    xend = data(end,2); % Extraire le x final
    zend = data(end,3); % Extraire le z final
    x_th = ((1.6022e-19*5/1.6726e-27)*4e5*cos(1.6022e-19*5*i/1.6726e-27)-(1.6726e-27*4e5/(1.6022e-19*5))); % TODO: Entrer la vraie solution analytique a tfin
    z_th = ((4e5*1.6726e-27/(1.6022e-19*5))*sin(1.6022e-19*5*i/1.6726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    error(i) = sqrt((xend-x_th)^2+(zend-z_th)^2); % erreur sur la position finale
end

figure
loglog(dt, error, 'k+')
xlabel('\Delta t')
ylabel('Maximum de l''erreur sur la position')
grid on

%Convergence de l'erreur sur la position

nsteps_num  = [100, 200, 400, 800, 1600, 3200, 6400];
xfin_ana = (-((37.0 - sqrt(111.0)) / (-34.0))*1e-6);
xfin_num = [7.75697638819292e-007,
            7.76455503400251e-007,
            7.76990057899759e-007,
            7.7736157050707e-007,
            7.77614331728026e-007,
            7.77780741836893e-007 ]
error_xfin = xfin_num-xfin_ana;
figure
plot(nsteps_num, abs(error_xfin),'k+-', "Color", "blue")
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('N_{steps}')
ylabel('Error on x_{fin}')
grid on

error = zeros(1,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    v_xend = data(end,2); % Extraire le x final
    v_zend = data(end,3); % Extraire le z final
    v_x_th = ((-1.6022e-19*1.6022e-19*5*5*4e5/(1.62726e-27*1.6726e-27))*sin(1.6022e-19*5*i/1.62726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    v_z_th = (4e5*cos(1.6022e-19*5*i/1.62726e-27)); % TODO: Entrer la vraie solution analytique a tfin
    error(i) = sqrt((v_xend-v_x_th)^2+(v_zend-v_z_th)^2); % erreur sur la position finale
end

figure
loglog(dt, error, 'k+')
xlabel('\Delta t')
ylabel('Maximum de l''erreur sur la vitesse')
grid on

%-----------------Fin de l'exo
%2.3.a.i--------------------------------------

