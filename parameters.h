//problem geometry, mesh control
#define DIMS 4 //2
#define FEOrder 2

#define problemHeight 1.0 //[1]
//#define problemWidth 1.0   //[2]
#define problemLength 1.0//[0]

#define globalRefinementFactor 0
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 0.005 //0.01
#define TotalTime 20004*TimeStep

//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 80 //70
#define YSubRf 80 //8
#define ZSubRf 1 //70

//Parameters for RBM
#define PRno 0.71
//#define GAMMA 6.0
#define RAno 10000.0
#define intT 0.05
//#define BIno 0.2
//#define MAno 92.0
//#define gravity 9.8
//#define Thot 0.135
//#define Tcold -0.115

/*
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
#define VV 0.0085

//Thermal Conductivity
#define KK 11.82+(1.06e-02)*T_conv[q]

#define KKS 60.0

//Heat Specific Capacity C
#define CC 330.9+(0.5653)*T_conv[q]-(4.015e-04)*T_conv[q]*T_conv[q]+(9.465e-08)*T_conv[q]*T_conv[q]*T_conv[q]

#define CCS 790.0

//Density of Material
#define RHO 7800.0 //4420.0

//Ambient Temperature
#define Tamb 301.3

//Heat Transfer coefficient
#define HH 10.0

//Latent Heat
#define LATENT 2.72e+05

//Radiation parameters
#define SIG (5.67037*pow(10.0,-8))
#define em 0.7  

//Laser parameters
#define DD 2.0
#define BB 2.0
#define spotRadius (0.6/1000.0)
#define PP 275.0
#define LAYER (0.254/1000.0)
#define ABSORB 0.99
*/




