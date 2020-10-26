% pour avoir une police pareil que sur Latex
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 16);

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

nsimul = 5; % Nombre de simulations a faire

nsteps = [785 1624 4642 6952 10000]; % Nombre d'iterations entier de 10^2 a 10^4  MODIFIER SELON VOS BESOINS
tfin = 1e-8; % TODO: Verifier que la valeur de tfin est la meme que dans le fichier input MODIFIER SELON VOS BESOINS
dt = tfin ./ nsteps;

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

lw = 2;

figure(1)
hold on
 for i=1:nsimul
    data = load(output{i});
    t = data(:,1);
    E = data(:,6);
    E0 = data(1,6);;
    plot(t,E-E0, 'linewidth', lw);
 end
    xlabel('t [s]');
    ylabel('$\Delta E$ [J]');
    grid on
