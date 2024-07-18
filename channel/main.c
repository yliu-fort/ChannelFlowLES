/**
## Flow over an obstacle

On this page we embed a free slip mountain into a domain that is filled with an viscous fluid. The free-slip condition is not free of rotation along a curved boundary. Therefore, simple potential flow theory does not apply and we solve the equations for fluid motion.  
*/
//#include <unistd.h>
//#include <stdlib.h>
//#include <stdio.h>

#include "navier-stokes/centered.h"
#if dimension == 3
#include "lambda2.h"
#endif

int kVerbose = 1;
#define MU (1.5e-5)
//#include "oystein/adapt_wavelet_leave_interface.h"
#include "sander/output_htg.h"

// Fast noise library https://github.com/Auburn/FastNoiseLite.git
#define FNL_IMPL
#include "FastNoiseLite.h"

#define MAX_LEVEL (8)
#define MIN_LEVEL (5)

// INPUT PARAMETERS
#define DENSITY (1.0)
#define CHANNEL_HEIGHT (1.0)
#define U_TAU (1.0)
#define RE_TAU (180.0)
#define UBAR (1.0)
#define DPDX (1.0)

#define Z0 (0.0001)

#define KAPPA (0.41)  // von Kármán constant
#define B (5.2)       // Log-law intercept (commonly used value for smooth walls)
scalar yp[], wallcorr[], utau[], roughness[];

/**
We set a free-slip boundary for the x and y components.
*/
scalar wall_indicator[];

u.n[top] = dirichlet (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
#if dimension > 2
u.r[bottom] = dirichlet (0);
#endif

wall_indicator[bottom] = dirichlet (1);

/**
Initialise pressure gradient force and viscosity vars.
*/
face vector muc[];
face vector av[];

int main(){
  L0 = 1.0;
  N = 1<<MIN_LEVEL;
  //kTurbConstants.k_0 = 0.001;
  //kTurbConstants.omega_0 = 420.0;
  a = av;
  mu = muc;

  CFL = 0.2;
  DT = 0.01;
  origin(-L0/2, 0, -L0/2);

  init_grid(N);
  periodic(right);
#if dimension > 2
  periodic(back);
#endif

  //TOLERANCE = 1e-5;
  run();
}

event properties(i=0){
  foreach_face(x)
    av.x[] = fm.x[]*DPDX;
  foreach_face()
    muc.x[] = fm.x[]*MU;

  boundary((scalar*){av, muc});
}

/**
We use particles to check if the flow does not penetrate the embedded boundary. 
*/
event init (t = 0){
/* The pressure gradient required to drive the flow is introduced as a
   source term */
   if (!restore ("restart")) {
  do{
    foreach(){
      //u.x[] = 1.842*UBAR*y*(2.0-y);
      u.x[] = UBAR;
      u.y[] = 0.0;
      u.z[] = 0.0;
      p[] = 0.;

      wall_indicator[] = 0.0;
    }

    boundary(all);
  }
#if TREE
    while (
    #if dimension == 2
    adapt_wavelet((scalar *){u,wall_indicator},(double[]){2e-3,2e-3, 0.001},MAX_LEVEL,MIN_LEVEL).nf
    #elif dimension == 3
    adapt_wavelet((scalar *){u,wall_indicator},(double[]){5e-3,5e-3,5e-3, 0.001},MAX_LEVEL,MIN_LEVEL).nf
    #endif
);
#else
    while (0);
#endif
   }
}

#if TREE
event adapt (i++) {
  #if dimension == 2
    adapt_wavelet((scalar *){u,wall_indicator},(double[]){2e-3,2e-3, 0.001},MAX_LEVEL,MIN_LEVEL);
  #elif dimension == 3
    adapt_wavelet((scalar *){u,wall_indicator},(double[]){5e-3,5e-3,5e-3, 0.001},MAX_LEVEL,MIN_LEVEL);
  #endif
}
#endif

double wall_function(double y_plus) {
    if (y_plus < 5.0) {
        return y_plus;
    }else if (5.0 <= y_plus <= 30.0){
        double gamma = (y_plus - 5.0) / 25.0;
        return pow((1.0 - gamma) * pow(5.0, 4.0) + gamma * pow((1.0 / KAPPA) * log(30.0) + B, 4.0), 0.25);
    }else if (y_plus > 30.0) {
        return (1.0 / KAPPA) * log(y_plus) + B;
    } else {
        return 0.0;  // Avoid log of zero or negative values
    }
}


// Function to return the fractional part of a number
double fract(double x) {
    return x - floor(x);
}

// https://www.shadertoy.com/view/4djSRW
// Noise function to generate a number in [0, 1]
double hash12(double p1, double p2) {
  double px = fract(p1 * .1031);
  double py = fract(p2 * .1031);
  double pz = fract(p1 * .1031);
  double q = px * (py + 33.33) + py * (pz + 33.33) + pz * (px + 33.33);

  px += q;
  py += q;
  pz += q;
  return fract((px + py) * pz);
}

double unoise(double p1, double p2) {
  return hash12(hash12(p1, 3.14159) * 27179.0, hash12(2.7179, p2) * 31415.9);
}


event wall_correction(i++){
  //fnl_state pnoise = fnlCreateState();
  //pnoise.noise_type = FNL_NOISE_PERLIN;
  //pnoise.fractal_type = FNL_FRACTAL_FBM;
  //pnoise.frequency = 30.0f;
  //pnoise.octaves = 5;

  foreach_boundary (bottom){
    double wall_dist = Delta/2.0; // bottom boundary is wall
    #if dimension == 2
      //double roughness_corr = fnlGetNoise2D(&pnoise, x, 0.0);
      double roughness_corr = unoise(x, 0.0);
    #elif dimension == 3
      //double roughness_corr = fnlGetNoise2D(&pnoise, x, z);
      double roughness_corr = unoise(x, z);
    #endif
    #if dimension == 2
      double actual_roughness = (1.0 + roughness_corr * (roughness_corr > 0.0 ? 1.0 : 0.8))*Z0;
      double u_bar = fabs(u.x[] - u.x[0,1]);
    #elif dimension == 3
      double actual_roughness = (1.0 + roughness_corr * (roughness_corr > 0.0 ? 1.0 : 0.8))*Z0;
      double u_bar = sqrt(sq(u.x[] - u.x[0,1]) + sq(u.z[] - u.z[0,1]));
    #endif
    double tau_w = -sq(KAPPA/log(wall_dist/actual_roughness))*sq(u_bar);
    double u_tau = sqrt(fabs(tau_w)/DENSITY);
    double y_plus = u_tau*wall_dist/(MU/DENSITY);

    yp[] = y_plus;
    utau[] = u_tau;
    roughness[] = actual_roughness;

    if (u_bar != 0) {  // Ensuring we don't divide by zero
        double correction_factor = min(1.0, wall_function(y_plus) * u_tau / u_bar);
        u.x[] *= correction_factor;
        #if dimension == 3
        u.z[] *= correction_factor;  // Assuming adjustment needed in both x and z directions
        #endif

        wallcorr[] = correction_factor;
    }
  }
}

/**
To help the accurate evaluation of the viscous term, we enforce the time step according to a limit cell-Peclet number. 
*/
event stability (i++) {
    CFL = 0.2;
}
/**
We generate a movie displaying the flow via $u_x$ and the tracers, these do not penetrate the boundary. 

![The particles follow the flow along the boundary](embed-freeslip/movie.mp4)
*/
event logfile_ke (i++) {
  //static FILE * fpke = fopen("ke.dat", "w");
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  
  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= MU/vol;

  if (i == 0)
    fprintf (stderr, "t dissipation energy Reynolds\n");
  fprintf (stderr, "%g %g %g %g\n",
	   t, vd, ke, 2./3.*ke/MU*sqrt(15.*MU/vd));
}

event logfile_mgp (i+=20) {
  if (i == 0)
    fprintf (stderr,
	     "i dt mgp.i mgp.resb mgp.resa mgp.sum mgp.nrelax TOL TOL/sq(dt)\n");
  fprintf (stderr, "%d %g %d %d %g %g %g %g %d %d %g %g \n",
	   i,dt, mgp.i,mgu.i,mgp.resb,mgu.resb,mgp.resa,mgu.resa,mgp.nrelax,mgu.nrelax,TOLERANCE,TOLERANCE/(dt*dt));
}

event snapshot(t += 0.05;t<=100.0) {
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  dump(file="restart", unbuffered=false);
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  fprintf (stderr, "%d : Dump restart file.\n", i);
}

event write2htg (t += 0.05;t<=100.0) {
  scalar vortex_core[];
  #if dimension == 2
    vorticity(u, vortex_core);
  #elif dimension == 3
    lambda2(u, vortex_core);
  #endif
  
  char fname[50];
  sprintf(fname, "result.%06d", i);
  scalar* output_scalars = {p, vortex_core, yp, wallcorr, utau, roughness};
  //rsm_model_output(&output_scalars);
  output_htg(output_scalars,(vector *){u}, ".", fname, i, t);

  fprintf (stderr, "write output to %s\n", fname);
}
