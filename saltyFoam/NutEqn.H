{
    alphat = turbulence->nut()/Prnt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Prn + alphan);

    fvScalarMatrix NuEqn
    (
        fvm::ddt(Nu)
      + fvm::div(phi, Nu)
      - fvm::laplacian(alphaEff, Nu)
     ==
       fvOptions(Nu)
    );

    NuEqn.relax();

    fvOptions.constrain(NuEqn);

    NuEqn.solve();

    radiation->correct();

    fvOptions.correct(Nu);

}
