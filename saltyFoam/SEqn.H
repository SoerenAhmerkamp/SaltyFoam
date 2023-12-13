{
    alphas = turbulence->nut()/Prst;
    alphas.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Prs + alphas); // turbulence->nu()/Pr + 

    fvScalarMatrix SEqn
    (
        fvm::ddt(S)
      + fvm::div(phi, S)
      - fvm::laplacian(alphaEff, S)
     ==
        fvOptions(S)
    );

    SEqn.relax();

    fvOptions.constrain(SEqn);

    SEqn.solve();

    radiation->correct();

    fvOptions.correct(S);

}
