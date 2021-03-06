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
nsteps = round(logspace(2,5,nsimul))

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
[x_th, z_th] = sol_anal_pos(1.6726e-27,1.6022e-19, 4e5, 5.0, 0, 7.e-8);
error_pos_max = zeros(1, nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    xend = data(end,2); % Extraire le x final
    zend = data(end,3); % Extraire le z final
    %x_th = 1.4768e13;
    %z_th = -1.9101e14;
    error_pos_max(i) = max(abs(xend-x_th), abs(zend-z_th)); % Maximum de l'erreur sur la position -> CORRECT
    
end
%Calcul de la pente de l'erreur sur la position,  on prend les dernières
%valeurs pour être plus précis
pente_pos = pente(dt(:, end-5:end), error_pos_max(:, end-5:end))

lw=1; fs=16;
figure('Name', [filename ': Convergence numérique de l''erreur sur la position'])
plot(dt, error_pos_max, 'b+-','linewidth',lw)
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Delta t')
ylabel('Convergence de l''erreur sur position')
grid on
 
 %x = x_th %Pour écrire dans la fenêtre (sans ;)
 %z = z_th

 

%Erreur max sur la vitesse
[v_xth, v_zth] = sol_anal_v(1.6726e-27,1.6022e-19, 4e5, 5.0, 0, 7.e-8);
error_vit_max = zeros(1, nsimul);

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    v_xend = data(end,4); % Extraire le v_x final
    v_zend = data(end,5); % Extraire le v_z final
    
    error_vit_max(i) = max(abs(v_xend-v_xth), abs(v_zend-v_zth)); % Maximum de l'erreur sur la position -> CORRECT
    
end
%Calcul de la pente de l'erreur sur la vitesse, on prend les dernières
%valeurs pour être plus précis
pente_vit = pente(dt(:,end-5:end), error_vit_max(:,end-5:end))

lw=1; fs=16; ms=6;
figure('Name', [filename ': Convergence numérique de l''erreur sur la vitesse'])
plot(dt, error_vit_max, 'b+-','linewidth',lw)
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Delta t')
ylabel('Convergence de l''erreur sur v')
grid on

%v_x = v_xth
%v_z = v_zth

%Erreur sur l'Energie mécanique 


delta_E_mec = zeros(1, nsimul);
[n] = n_steps_limit(1.6726e-27,1.6022e-19, 5.0, 1.e-8) %Valeur limite de n_steps pour laquelle le schéma d'EC est stable
[E_mec_0] = E_mec_i(1.6726e-27, 4e5)
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    dt(i) = data(2,1)-data(1,1); % Delta t = t_2-t_1
    E_mec_t = data(end,6); % Extraire E_mec(t)
    delta_E_mec(i) = E_mec_t - E_mec_0; % Delta E_mec
end

figure('Name', [filename ': Convergence numérique de l''erreur sur l''énergie mécanique'])
plot(dt,delta_E_mec, '+', 'linewidth',lw)
set(gca,'fontsize',fs)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\Delta t')
ylabel('\Delta E_{mec}')
grid on


the = x_th
worst = z_th
idea = v_xth
ever = v_zth

%Fonction qui calcule la position théorique
function [x,z] = sol_anal_pos(m, q, v_0, B_0, E_0, t_fin)
   
    
    %x_th selon x
    A = (q*B_0)./m;
    C1 = v_0./A;
    C2 = E_0/(A.*B_0);
    
    
     x = C1.*cos(A.*t_fin)+ C2.*sin(t_fin*A) - (E_0.*t_fin./B_0) - (v_0./A);
   
     %z_th selon z
     
     z =  -C2.*cos(t_fin*A)+ C1.*sin(t_fin*A) + (E_0/(A.*B_0));
   

end


%Fonction qui calcule la vitesse v_th
function [v_x,v_z] = sol_anal_v(m, q, v_0, B_0, E_0, t_fin)
    %v_th selon x
    A = (q.*B_0)./m;
     v_x = -v_0.*sin(A*t_fin) + (E_0/B_0).*cos(A*t_fin) - (E_0/B_0);
   
     %v_th selon z 
     
     v_z = (E_0/B_0).*sin(A*t_fin) + v_0.*cos(A*t_fin);
   

end

%Fonction qui calcule le nsteps limite pour Euler Cromer
function [n] = n_steps_limit(m, q, B_0, t_fin) 
   w = q.*B_0./m;
   n = w.*t_fin./2.0;

end

%Fonction qui calcule l'E_mec initiale
function [E_mec_0] = E_mec_i(m, v_0) 
  E_mec_0 = 1.0./2.0 *m*v_0*v_0;
end

function [pent] = pente(nsteps, error)
    pent = log( flip(error,2) ./ error) ./ log(flip(nsteps,2) ./ nsteps);
end

