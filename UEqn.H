// Solve the Momentum equation

fvVectorMatrix UEqn
(
    fvm::ddt(U) 
	+ fvm::div(phi, U)
	+ turbulence->divDevReff(U)
	- rho*E/rho0
	==
    fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
