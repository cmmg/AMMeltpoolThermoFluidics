//problem geometry, mesh control
#define DIMS 4 //2
#define FEOrder 2

#define problemWidth 25.0  //[2]
#define problemHeight 4.1   //[1]
#define problem_Length 0.01 //[0]

#define globalRefinementFactor 3
#define maxRefinementLevel (globalRefinementFactor+2)
#define minRefinementLevel (globalRefinementFactor)

//time step controls
#define TimeStep 5.0e-4
#define TotalTime 5000*TimeStep

//output controls
#define outputFileName "solution"

//subdivisons
#define XSubRf 1
#define YSubRf 1
#define ZSubRf 1

//Kinematic viscosity
#define nu 1.0



