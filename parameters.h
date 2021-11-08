//problem geometry, mesh control
#define DIMS 6 //2
#define FEOrder 2

#define problemWidth  25.0 //[2]
#define problemHeight 25.0   //[1]
#define problemLength 150.0   //[0]

#define globalRefinementFactor 1
#define maxRefinementLevel (globalRefinementFactor+3)
#define minRefinementLevel (globalRefinementFactor-1)

//time step controls
#define TimeStep 1.0e-00
#define TotalTime 217*TimeStep
#define TOLERANCE 1.0e-08


//output controls
#define outputFileName "solution"
#define PSTEPS 1

//subdivisons
#define XSubRf 6 //53
#define YSubRf 1 //106
#define ZSubRf 1

//Material parameters of ss316
//Non-dimensional numbers
#define PRno (9.40e-02)
#define PEno (6.28e-01)
#define GRno (7.32e-09)
#define MAno (-1.71e+02)
#define STEno (7.21e-02)
#define TCONST (2.79e-02)
#define BIno (5.48e-06)
#define TSBOno (1.45e-04)
#define DAno (1.39e-03)
#define rs (5.0)
#define th (1.0)
#define Vs (1.0)
#define Qlaser (1.07e+02)
#define Tinit (3.62e-02)
#define emsvty (0.9)
#define DD (2.0)
#define TSS (9.72e-1)
#define Tamb (2.10e-01)
