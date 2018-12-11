%Partie 3: Méthodes de gradient accélérés
%Fonction de l'erreur
function E = erreur(A,t,y)
a = A(1);
b = A(2);
yh = a*(1-exp(b*t));
E = log(sum((y-yh).ˆ2));

%Gradient de l'erreur
function gradE = grad erreur(A,t,y)
a = A(1);
b = A(2);
yh = a*(1-exp(b*t)); % Vecteur des valeur estimées de y
dya = (1-exp(b*t)); % Dérivée de yh par rapport à a
dyb = a*(-t.*exp(b*t)); % Dérivée de yh par rapport à b
S = sum((y-yh).ˆ2); % Sommes des carrés des erreurs
daS = 2*(yh-y)’*dya; % Dérivée de S par rapport à a
dbS = 2*(yh-y)’*dyb; % Dérivée de S par rapport à b
gradE = [daS ; dbS]/S; % Gradient de l’erreur logarithmique

%Visualisation de la fonction sous forme de surface et de contours
load partie3
a = 0.5:0.01:2;
b = -2:0.01:0.5;
[A,B]=meshgrid(a,b);
[n,m]=size(A);
for i=1:n
  for j=1:m
    aa = [A(i,j) ; B(i,j)];
    ERR(i,j)=erreur(aa,t,y);
  end
end
mesh(A,B,ERR);
figure
contour(a,b,ERR,100);
xlabel(’x’);ylabel(’y’)
%Les courbes de niveaux ne sont pas convexes donc la fonction n'est pas convexe globalement.
%Elle est convexe seulement localement au voisinage du minimum.

%Méthode du gradient à pas fixe
function [X,i] = pas constant(A0,t,y,rho,eps)
X = A0;
g = grad erreur(A0,t,y);
A1 = A0 - rho*g;
X = [X A1];
i = 0;
while norm(A1-A0) > eps && i<50000
  i=i+1;
  A0 = A1;
  g = grad erreur(A0,t,y);
  A1 = A0 - rho*g;
  X = [X A1];
end

%Méthode d'inertie "Méthode du moment"
function [X,i] = moment(A0,v0,t,y,rho,eta,eps)
X = A0;
g = grad erreur(A0,t,y);
v1 = eta*v0+rho*g;
A1 = A0 - v1;
X = [X A1];
i = 0;
while norm(A1-A0) > eps && i<50000
  i=i+1;
  A0 = A1;
  v0 = v1;
  g = grad erreur(A0,t,y);
  v1 = eta*v0+rho*g;
  A1 = A0 - v1;
  X = [X A1];
end

%Méthode du gradient accéléré de Nesterov
function [X,i] = nesterov(A0,v0,t,y,rho,eta,eps)
X = A0;
g = grad erreur(A0-eta*v0,t,y);
v1 = eta*v0+rho*g;
A1 = A0 - v1;
X = [X A1];
i = 0;
while norm(A1-A0) > eps && i<50000
  i=i+1;
  A0 = A1;
  v0 = v1;
  g = grad erreur(A0-eta*v0,t,y);
  v1 = eta*v0+rho*g;
  A1 = A0 - v1;
  X = [X A1];
end

%Visualisation des trajectoires
load partie3
eps = 1e-6;
rho = 0.00035;
eta = 0.09;
A0 = [1;0];
v0 = grad erreur(A0,t,y);
X = pas constant(A0,t,y,rho,eps);
X1 = moment(A0,v0,t,y,rho,eta,eps);
X2 = nesterov(A0,v0,t,y,rho,eta,eps);
open(’partie3-contours-1.fig’)
hold on plot(X(1,:),X(2,:));
plot(X1(1,:),X1(2,:))
plot(X2(1,:),X2(2,:),’k’);
legend(’Courbes de niveau’,’Pas constant’,’Inertie’,’Gradient accéléré’);

%Nombre d'itérations en fonction de -log(epsilon)
load partie3
e = [3:9];
eps = 10.ˆ(-e);
rho = 0.00035;
eta = 0.09;
A0 = [1;0];
v0 = grad erreur(A0,t,y);
for k=1:length(eps)
  [X,i] = pas constant(A0,t,y,rho,eps(k));
  N(k)=i;
  [X1,i1] = moment(A0,v0,t,y,rho,eta,eps(k));
  N1(k)=i1;
  [X2,i2] = nesterov(A0,v0,t,y,rho,eta,eps(k));
  N2(k)=i2;
end
plot(e,[N’ N1’ N2’]);
grid
legend(’Pas constant’,’Inertie’,’Gradient accéléré’);
xlabel(’-log(nepsilon)’);
ylabel(’Nombre d’’itérations’);
%Les méthodes d'inertie et de gradient accéléré ont à peu près la même performance et
%sont plus performantes que la méthode du gradient à pas fixe

%Nombre d'itérations et le temps de calcul pour chaque méthode
load partie3
eps = 1e-6;
rho = 0.00035;
eta = 0.09;
A0 = [-2;-1];
v0 = grad erreur(A0,t,y);
tic
[X,i] = pas constant(A0,t,y,rho,eps);
display([’Gradient `a pas fixe j Nombre d’’itérations : ’ num2str(i) ’ j Temps
d’’exécution:’ num2str(toc) ]);
tic
[X1,i1] = moment(A0,v0,t,y,rho,eta,eps);
display([’Méthode d’’inertie j Nombre d’’itérations : ’ num2str(i1) ’ j Temps
d’’exécution:’ num2str(toc) ]);
tic [X2,i2] = nesterov(A0,v0,t,y,rho,eta,eps);
display([’Gradient accéléré de Nesterov j Nombre d’’itérations : ’ num2str(i1) ’ j
Temps d’’exécution:’ num2str(toc)]);s
%La surface au voisinage de la solution contient un creux. Le gradient à pas fixe présente
%des oscillations car le pas choisi ne permet pas de s'adapter au creux.

%Animation
load partie3
eps = 1e-6;
rho = 0.00035;
eta = 0.09;
A0 = [0.1;-0.2];
v0 = grad erreur(A0,t,y);
[A,i2] = nesterov(A0,v0,t,y,rho,eta,eps);
n = length(A);
for i=10:10:n
  plot(t,y,’.’);
  hold on
  a = A(1,i);
  b = A(2,i);
  yh = a*(1-exp(b*t));
  plot(t,yh,’r’);
  hold off
  xlabel(’Temps’);ylabel(’y’);
  legend(’Valeurs mesur´ees’,’Valeurs estim´ees’,’location’,’SE’);
  % ’location’,’SE’ permet de placer la l´egende en bas `a droite de la figure
  title([’It´eration : ’ num2str(i)]);
  set(gcf,’color’,’w’); % Pour avoir un fond blanc
  getframe();
end
