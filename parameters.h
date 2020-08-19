//problem geometry, mesh control
#define DIMS 5 //2
#define FEOrder 2

#define problemWidth 5.0e-03  //[2]
#define problemHeight 5.0e-03   //[1]
#define problemLength 20.0e-03 //[0]

#define globalRefinementFactor 3
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 1.0e-3
#define TotalTime 20000*TimeStep

//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 1
#define YSubRf 1
#define ZSubRf 1



//Material parameters

//viscosity
#define mu 7.0e-03

//surface tension grad
#define dGammadT -0.4e-03

//expansion coeff
#define BETA 5.85e-05

//PDAS in micron
#define PDAS 0.5 


//Liquidus Temperature
#define TLL 1928.0

//Solidus Temperature
#define TSS 1878.0

//scan speed
#define VV 0.0033*10.0

//Thermal Conductivity
#define KK 34.0

//Heat Specific Capacity C
#define CC 830.0

//Density of Material
#define RHO 4420.0

//kinematic viscosity
#define nu mu/RHO

//Ambient Temperature
#define Tamb 300.0

//Heat Transfer coefficient
#define HH 24.0

//Latent Heat
#define LATENT (2.0*pow(10.0,5))

//Radiation parameters
#define SIG (5.67037*pow(10.0,-8))
#define em 0.9  

//Laser parameters
#define DD 0.5
#define BB 2.0
#define spotRadius (0.5/1000.0)
#define PP 400.0
#define LAYER (0.1/1000.0)
#define ABSORB 1.0


//gravity
#define gravity 9.8
