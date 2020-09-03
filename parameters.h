//problem geometry, mesh control
#define DIMS 5 //2
#define FEOrder 2

#define problemWidth  //[2]
#define problemHeight 6.35/100.0   //[1]
#define problemLength 3.175/100.0 //[0]

#define globalRefinementFactor 1
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 0.2
#define TotalTime 500*TimeStep

//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 53
#define YSubRf 106
#define ZSubRf 1

//Material parameters

//viscosity
#define mu 1.81e-03

//surface tension grad
#define dGammadT -3.7e-04

//expansion coeff
#define BETA 1.2e-04

//PDAS in micron
#define PDAS 0.5

//Liquidus Temperature
#define TLL 305.78 //1928.0

//Solidus Temperature
#define TSS 302.78 //1878.0

//scan speed
#define VV 16.7e-04

//Thermal Conductivity
#define KK 32.0

#define KKL 32.0

//Heat Specific Capacity C
#define CC 381.5

#define CCL 381.5

//Density of Material
#define RHO 6093.0 //4420.0

#define RHOL 6093.0 //4420.0

//kinematic viscosity
#define nu mu/RHOL

//kinematic viscosity
#define mus 1.0e+04

//Ambient Temperature
#define Tamb 301.3

//Heat Transfer coefficient
#define HH 100.0

//Latent Heat
#define LATENT 80160.0 //(2.72*pow(10.0,5))

//Radiation parameters
#define SIG (5.67037*pow(10.0,-8))
#define em 0.9  

//Laser parameters
#define DD 2.0
#define BB 2.0
#define spotRadius (1.5/1000.0)
#define PP 4000.0
#define LAYER (0.1/1000.0)
#define ABSORB 1.0


//gravity
#define gravity 10.0


