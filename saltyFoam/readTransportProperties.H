singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);
dimensionedScalar Prs("Prs", dimless, laminarTransport);
dimensionedScalar Prst("Prst", dimless, laminarTransport);
dimensionedScalar Prnt("Prnt", dimless, laminarTransport);
dimensionedScalar Prn("Prn", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);
