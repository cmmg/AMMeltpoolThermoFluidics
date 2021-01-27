//problem geometry, mesh control
#define DIMS 6 //2
#define FEOrder 2

#define problemHeight 1.61e-04 //[1]
#define problemWidth 2.0e-04   //[2]
#define problemLength 6.0e-04//[0]

#define globalRefinementFactor 0
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 1.0e-06
#define TotalTime 554*TimeStep
#define TOLERANCE 1.0e-08

//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 60 //60
#define YSubRf 20 //15
#define ZSubRf 20 //35

//Material parameters of ss316

//viscosity
#define mu 3.0e-03

//surface tension grad
#define dGammadT -2.13e-04

//expansion coeff
#define BETA 9.541e-05

//PDAS in micron
#define PDAS 0.5

//Liquidus Temperature
#define TLL 868.0

//Solidus Temperature
#define TSS 743.0

//melting range
#define deltaT 125.0

//scan speed
#define VV 1.0

//Thermal Conductivity
#define KKS 72.0

#define KKL 82.9

//Heat Specific Capacity C
#define CCL 1230.0

#define CCS 1014.0

//Density of Material
#define RHOS 1809.5//4420.0

//Density of Material
#define RHOL (1731.2 -0.2208*(T_conv[q]-TLL)) //4420.0


//Ambient Temperature
#define Tamb 293.15

//Heat Transfer coefficient
#define HH 10.0 //10.0

//Latent Heat
#define LATENT 3.73e+05

//Radiation parameters
#define SIG (5.67037*pow(10.0,-8))
#define em 0.18  

//Laser parameters
#define DD 2.0
#define BB 2.0
#define spotRadius (0.08/1000.0)
#define PP 100.0
#define LAYER (0.04/1000.0)
#define ABSORB 0.18
#define porosity 0.475

//gravity
#define gravity 9.8


