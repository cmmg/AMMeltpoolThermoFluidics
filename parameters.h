//problem geometry, mesh control
#define DIMS 3 //2
#define FEOrder 2

#define problem_Width 0.003  //[2]
#define problem_Height 0.003   //[1]
#define problem_Length 0.01 //[0]

#define globalRefinementFactor 1
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 0.01
#define TotalTime 1500*TimeStep

//output controls
#define outputFileName "solution"

//Heat Diffusion parameters

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
#define LAYER (0.5/1000.0)
#define ABSORB 0.99

//subdivisons
#define XSubRf 8
#define YSubRf 2
#define ZSubRf 2


