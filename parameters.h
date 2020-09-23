//problem geometry, mesh control
#define DIMS 5 //2
#define FEOrder 2

#define problemWidth  (5.0/1000.0) //[2]
#define problemHeight (5.0/1000.0)   //[1]
#define problemLength (20.0/1000.0)//[0]

#define globalRefinementFactor 0
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 0.1
#define TotalTime 250*TimeStep

//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 200 //53
#define YSubRf 50 //106
#define ZSubRf 1

//Material parameters of ss316

//viscosity
#define mu 7.0e-03

//surface tension grad
#define dGammadT -0.4e-03

//expansion coeff
#define BETA 5.85e-05

//PDAS in micron
#define PDAS 0.5

//Liquidus Temperature
#define TLL 1733.0

//Solidus Temperature
#define TSS 1693.0

//melting range
#define deltaT 40.0

//scan speed
#define VV 1.0e-03

//Thermal Conductivity
#define KK 11.82+(1.06e-02)*T_conv[q]

//#define KK 32.0

//Heat Specific Capacity C
#define CC 330.9+(0.5653)*T_conv[q]-(4.015e-04)*T_conv[q]*T_conv[q]+(9.465e-08)*T_conv[q]*T_conv[q]*T_conv[q]

//#define CC 381.5

//Density of Material
#define RHO 7800.0 //4420.0

//Ambient Temperature
#define Tamb 301.3

//Heat Transfer coefficient
#define HH 100.0

//Latent Heat
#define LATENT 2.72e+05

//Radiation parameters
#define SIG (5.67037*pow(10.0,-8))
#define em 0.9  

//Laser parameters
#define DD 2.0
#define BB 2.0
#define spotRadius (1.5/1000.0)
#define PP 2.0*90.0
#define LAYER (0.1/1000.0)
#define ABSORB 1.0


//gravity
#define gravity 9.8


