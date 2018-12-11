%Partie 2: Optimisation avec contrainte : la régression linéaire multiple
%Calcul du gradient
function e = gradE(a,lambda,A,b)
e = 2*A*a-2*b+2*lambda*a;
%Calcul de la hessienne
function H = hessE(a,lambda,A)
H = 2*A+2*lambda*eye(length(a));
%Méthode de Newton
function [a1,i] = newton(a0,lambda,A,b,eps)
g = gradE(a0,lambda,A,b);
H = hessE(a0,lambda,A);
a1 = a0-Hng;
i=0;
while (norm(a0-a1)>eps)
  i=i+1;
  a0 = a1;
  g = gradE(a0,lambda,A,b);
  H = hessE(a0,lambda,A);
  a1 = a0-Hng;
end

%Traçage de l'évolution des coefficients en fonction de lambda
load partie2
a0 = zeros(7,1);
A = X’*X;
b = X’*y;
a coef = [];
esp = 1e-3;
L = [1/1000:1/1000:1 2:500]; % de 0.0001 à 1 par pas de 0.0001 puis de 2 à 500 par pas de 1
for lambda = L
  a = newton(a0,lambda,A,b,eps);
  a coef = [a coef a];
end
plot(L,a coef’);
legend(’a 0’,’a 1’,’a 2’,’a 3’,’a 4’,’a 5’,’a 6’)
%La pénalisation est efficace car elle permet bien d'éviter le phénomèene d'explosion des coefficients.

%Traçage de l'évolution des variables y et y^opt(lambda) pour différentes valeurs de lambda :
plot(y,’r.’,’markersize’,10);
n = size(a coef,2);
hold on
for i=1:n
a = a coef(:,i);
h = X*a;
plot(yh);
end
axis([0 115 -9 10])
legend(’y (donn´ees)’,’y foptg(nlambda=1)’,...
’y foptg(nlambda=100)’,’y foptg(nlambda=200)’,...
’y foptg(nlambda=300)’,’y foptg(nlambda=400)’,’y foptg(nlambda=500)’)
%Plus lambda est grand, plus les courbes des valeurs estimées y^opt(lambda) s'éloignent de y.
