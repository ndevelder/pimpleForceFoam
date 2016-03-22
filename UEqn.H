// Solve the Momentum equation

// Add forcing term
const dimensionedScalar forcingTerm(transportProperties.lookup("forcingTerm"));
Info<< "Forcing Term: " <<forcingTerm <<endl;
dimensionedVector forcing("forcing", dimensionSet(0,1,-2,0,0,0,0),vector(forcingTerm.value(),0.0,0.0));

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  + forcing
);

Info<< "Added Forcing to NS Equation" <<endl;

if (oCorr == nOuterCorr - 1)
{
    if (mesh.solutionDict().relax("UFinal"))
    {
        UEqn().relax(mesh.solutionDict().relaxationFactor("UFinal"));
    }
    else
    {
        UEqn().relax(1);
    }
}
else
{
    UEqn().relax();
}

volScalarField rUA = 1.0/UEqn().A();

if (momentumPredictor)
{
    if (oCorr == nOuterCorr-1)
    {
        solve(UEqn() == -fvc::grad(p), mesh.solutionDict().solver("UFinal"));
    }
    else
    {
        solve(UEqn() == -fvc::grad(p));
    }
}
else
{
    U = rUA*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}
