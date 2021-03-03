//problem geometry, mesh control
#define DIMS 6 //2
#define FEOrder 2

#define problemWidth  0.5e-03 //[2]
#define problemHeight 0.5e-03   //[1]
#define problemLength 3.0e-03   //[0]

#define globalRefinementFactor 1
#define maxRefinementLevel (globalRefinementFactor+4)
#define minRefinementLevel (globalRefinementFactor-1)

//time step controls
#define TimeStep 1.0e-5
#define TotalTime 704*TimeStep
#define TOLERANCE 1.0e-08


//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 6 //53
#define YSubRf 1 //106
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
#define VV 0.575

//Thermal Conductivity
#define KKL 60.0 //11.82+(1.06e-02)*T_conv[q]

#define KKS 11.82+(1.06e-02)*T_conv[q]


//#define KK 32.0

//Heat Specific Capacity C
#define CCL 790.0 //330.9+(0.5653)*T_conv[q]-(4.015e-04)*T_conv[q]*T_conv[q]+(9.465e-08)*T_conv[q]*T_conv[q]*T_conv[q]

#define CCS 330.9+(0.5653)*T_conv[q]-(4.015e-04)*T_conv[q]*T_conv[q]+(9.465e-08)*T_conv[q]*T_conv[q]*T_conv[q]

//Density of Material
#define RHO 7800.0 //4420.0

//Ambient Temperature
#define Tamb 300.0

//Heat Transfer coefficient
#define HH 10.0

//Latent Heat
#define LATENT 2.72e+05

//Radiation parameters
#define SIG (5.67037*pow(10.0,-8))
#define em 0.9

//Laser parameters
#define DD 2.0
#define BB 2.0
#define spotRadius (0.1/1000.0)
#define PP 90.0
#define LAYER (0.02/1000.0)
#define ABSORB 0.9


//gravity
#define gravity 9.8


