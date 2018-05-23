function [u, b, C]=  lse_bfe_3Phase(u, img, b, Ksigma, KONE, nu, timestep, mu, epsilon)

% This code implements the level set method by Chuming Li (2011, IEEE
% transaction), by Pengwei Wu, guided by Tianye Niu.
% The code implemented is now three-phase (which can be extended to four or
% more phases).


N_class = 3;  % three-phase
KB1 = conv2(b, Ksigma, 'same');
KB2 = conv2(b.^2, Ksigma, 'same');
H1 = Heaviside(u(:,:,1), epsilon);
H2 = Heaviside(u(:,:,2), epsilon);
M(:,:,1) = H1 .* H2; % membership function 1
M(:,:,2) = H1 .* (1 - H2); % membership function 2
M(:,:,3) = (1 - H1); % membership function 3
C = updateC(img, KB1, KB2, M);

KONE_img = img .^ 2 .* KONE;
u = updateLSF(img, u, C, N_class, KONE_img, KB1, KB2, mu, nu, timestep, epsilon);

b = updateB(img, C, M,  Ksigma);
end


function u = updateLSF(img,u, C, N_class, KONE_img, KB1, KB2, mu, nu, timestep, epsilon)
u(:,:,1) = NeumannBoundCond(u(:,:,1));
Curv(:,:,1) = curvature_central(u(:,:,1));
H1 = Heaviside(u(:,:,1), epsilon);
Delta(:,:,1) = Dirac(u(:,:,1), epsilon);

u(:,:,2) = NeumannBoundCond(u(:,:,2));
Curv(:,:,2)=curvature_central(u(:,:,2));
H2 =  Heaviside(u(:,:,2), epsilon);
Delta(:,:,2) = Dirac(u(:,:,2), epsilon);

e = zeros([size(img),N_class]);
for kk = 1:N_class
    e(:,:,kk) = KONE_img - 2 * img .* C(kk) .* KB1 + C(kk) ^ 2 * KB2;
end

A1 = - Delta(:,:,1) .* (e(:,:,1) .* H2 + e(:,:,2) .* (1-H2) - e(:,:,3));
P1 = mu * (4*del2(u(:,:,1)) - Curv(:,:,1));
L1 = nu .* Delta(:,:,1) .* Curv(:,:,1);
u(:,:,1) = u(:,:,1) + timestep * (L1+P1+A1);   % update u1


A2 = - Delta(:,:,2) .* H1 .* (e(:,:,1) - e(:,:,2));
P2 = mu * (4 * del2(u(:,:,2)) - Curv(:,:,2));
L2=nu.*Delta(:,:,2) .* Curv(:,:,2);
u(:,:,2) = u(:,:,2) + timestep * (L2 + P2 + A2);   % update u2
end

% update A (u1 and u2) is done in the function above

function C=updateC(img, Kb1, Kb2, M)

N_class = size(M, 3);
for kk=1:N_class
    N2 = Kb1.*img.*M(:,:,kk);
    D2 = Kb2.*M(:,:,kk);
    sN2 = sum(N2(:));
    sD2 = sum(D2(:));
    C(kk) = sN2 / (sD2 + (sD2 == 0)); %#ok<AGROW>
end
end

function b = updateB(img, C, M, Ksigma)
PC1 = zeros(size(img));
PC2 = PC1;
N_class = size(M,3);
for kk = 1:N_class
    PC1 = PC1 + C(kk) * M(:,:,kk);
    PC2 = PC2 + C(kk) ^ 2 * M(:,:,kk);
end
KNm1 = conv2(PC1.*img, Ksigma, 'same');
KDn1 = conv2(PC2, Ksigma, 'same');

b = KNm1 ./ KDn1;
end


function f = Dirac(x, epsilon)    % function (12)
f = (epsilon / pi) ./ (epsilon ^ 2. + x .^ 2); % 2. for double conversion
end


function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow], [1 ncol]) = g([3 nrow-2], [3 ncol-2]);
g([1 nrow], 2:end-1) = g([3 nrow-2], 2:end-1);
g(2:end-1, [1 ncol]) = g(2:end-1, [3 ncol-2]);
end