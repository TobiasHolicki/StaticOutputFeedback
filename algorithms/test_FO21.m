function test_FO21()

[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC3');

sys.A=A;
sys.Bw=B1;
sys.Bu=B;
sys.Cy=C;
sys.Cz=C1;
sys.Dzu=D12;
sys.Dzw=D11;
sys.Dyw=D21;

opt.Kini = [  0.1250    0.1250   -0.1250    0.1250
   -0.1250   -0.1250    0.1250   -0.1250];
   
opt.Kini = [   -0.01133     -0.125      0.125      0.125; ...
    -0.006417     -0.125      0.125      0.125];

opt.echo=1;
opt.maxIt = 50;
opt.nc=0;
out = dof_hinf_FO21_c(sys,opt)

if out.feas
    out.L
end