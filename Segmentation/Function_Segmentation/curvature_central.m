function K_curvature = curvature_central(u)
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-20);
Nx = ux./normDu;
Ny = uy./normDu;
[nxx,~] = gradient(Nx);
[~,nyy] = gradient(Ny);
K_curvature = nxx+nyy;
end