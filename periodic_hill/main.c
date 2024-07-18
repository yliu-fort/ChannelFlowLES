/**
## Flow over an obstacle

On this page we embed a free slip mountain into a domain that is filled with an viscous fluid. The free-slip condition is not free of rotation along a curved boundary. Therefore, simple potential flow theory does not apply and we solve the equations for fluid motion.  
*/
#include "embed.h"
#include "navier-stokes/centered.h"

int kVerbose = 1;
#define MU (1e-5)
#include "oystein/adapt_wavelet_leave_interface.h"
#include "sander/output_htg.h"

#define MAX_LEVEL (8)
#define MIN_LEVEL (5)

// INPUT PARAMETERS
#define DENSITY (1.0)
#define CHANNEL_HEIGHT (1.0)
#define U_TAU (1.0)
#define RE_TAU (180.0)
#define UBAR (1.842)
#define DPDX (1.0)

/**
We set a free-slip boundary for the x and y components.
*/
u.n[top] = dirichlet (0);

u.n[embed] = dirichlet (0);
u.t[embed] = dirichlet (0);
#if dimension > 2
u.r[embed] = dirichlet (0);
#endif

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
#if dimension > 2
u.r[bottom] = dirichlet (0);
#endif

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

  //TOLERANCE = 1e-8;
  run();
}

event properties(i=0){
  foreach_face(x)
    av.x[] = fm.x[]*DPDX;
  foreach_face()
    muc.x[] = fm.x[]*MU;

  boundary((scalar*){av, muc});
}

double humpx(double xc)
{
    double x = 1000.0*xc/L0;
    double h = 0;

    //if (x > 198) x = 252 - x;

    if (x >= 0 && x < 9)
    {
        h =
            28
          + 6.775070969851E-03*x*x
          - 2.124527775800E-03*x*x*x;
    }
    else if (x >= 9 && x < 14)
    {
        h =
            25.07355893131
          + 0.9754803562315*x
          - 1.016116352781E-01*x*x
          + 1.889794677828E-03*x*x*x;
    }
    else if (x >= 14 && x < 20)
    {
        h =
            2.579601052357E+01
          + 8.206693007457E-01*x
          - 9.055370274339E-02*x*x
          + 1.626510569859E-03*x*x*x;
    }
    else if (x >= 20 && x < 30)
    {
        h =
            4.046435022819E+01
          - 1.379581654948E+00*x
          + 1.945884504128E-02*x*x
          - 2.070318932190E-04*x*x*x;
    }
    else if (x >= 30 && x < 40)
    {
        h =
            1.792461334664E+01
          + 8.743920332081E-01*x
          - 5.567361123058E-02*x*x
          + 6.277731764683E-04*x*x*x;
    }
    else if (x >= 40 && x < 54)
    {
        h =
            max
            (
                0,
                5.639011190988E+01
              - 2.010520359035E+00*x
              + 1.644919857549E-02*x*x
              + 2.674976141766E-05*x*x*x
            );
    }

    return h/1000.0;
}

double hump(double xc, double zc)
{
  double r0 = sqrt((xc+L0/2)*(xc+L0/2) + (zc+L0/2)*(zc+L0/2));
  double r1 = sqrt((xc+L0/2)*(xc+L0/2) + (L0/2-zc)*(L0/2-zc));
  double r2 = sqrt((L0/2-xc)*(L0/2-xc) + (zc+L0/2)*(zc+L0/2));
  double r3 = sqrt((L0/2-xc)*(L0/2-xc) + (L0/2-zc)*(L0/2-zc));
  double r = min(min(r0, r1), min(r2, r3));
    return humpx(r)-L0/(1<<(1+MAX_LEVEL));
}

/**
We use particles to check if the flow does not penetrate the embedded boundary. 
*/
event init (t = 0){
/* The pressure gradient required to drive the flow is introduced as a
   source term */

  do{
    scalar phi[];
    foreach_vertex(){
      #if dimension == 2
        phi[] = y - hump(x, L0/2);
      #elif dimension == 3
        phi[] = y - hump(x, z);
      #endif
    }
    fractions(phi, cs, fs);

    foreach(){
      u.x[] = UBAR*cs[];
      u.y[] = 0.0;
      u.z[] = 0.0;
      p[] = 0.;
    }

    boundary(all);
  }
#if TREE
    while (
    #if dimension == 2
    adapt_wavelet_leave_interface((scalar *){u},{cs},(double[]){2e-3,2e-3},MAX_LEVEL,MAX_LEVEL,MIN_LEVEL,1).nf
    #elif dimension == 3
    adapt_wavelet_leave_interface((scalar *){u},{cs},(double[]){5e-3,5e-3,5e-3},MAX_LEVEL,MAX_LEVEL,MIN_LEVEL,1).nf
    #endif
);
#else
    while (0);
#endif

  output_ppm(cs, n = 512, file = "cs.png", min = 0, max = 1);
}

#if TREE
event adapt (i++) {
  #if dimension == 2
    adapt_wavelet_leave_interface((scalar *){u},{cs},(double[]){2e-3,2e-3},MAX_LEVEL,MAX_LEVEL,MIN_LEVEL,1);
  #elif dimension == 3
    adapt_wavelet_leave_interface((scalar *){u},{cs},(double[]){5e-3,5e-3,5e-3},MAX_LEVEL,MAX_LEVEL,MIN_LEVEL,1);
  #endif
}
#endif


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
event logfile (i++) {
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
    fprintf (ferr,
	     "i dt mgp.i mgp.resb mgp.resa mgp.sum mgp.nrelax TOL TOL/sq(dt)\n");
  fprintf (ferr, "%d %g %d %d %g %g %g %g %d %d %g %g \n",
	   i,dt, mgp.i,mgu.i,mgp.resb,mgu.resb,mgp.resa,mgu.resa,mgp.nrelax,mgu.nrelax,TOLERANCE,TOLERANCE/(dt*dt));
}

event snapshot(t += 0.05;t<=5.0) {
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  dump(file="restart", unbuffered=false);
#if _MPI  
  MPI_Barrier(MPI_COMM_WORLD);  
#endif
  fprintf (stderr, "%d : Dump restart file.\n", i);
}

event write2htg (t += 0.05;t<=5.0) {
  scalar vort[];
  vorticity (u, vort);

  vector U[]; // to slove the compatibility issue of restore/dump
  foreach()
    foreach_dimension()
      U.x[] = u.x[];
  boundary((scalar *){U});
  
  char fname[50];
  sprintf(fname, "result.%06d", i);
  scalar* output_scalars = {cs, p, vort};
  //rsm_model_output(&output_scalars);
  output_htg(output_scalars,(vector *){U}, ".", fname, i, t);

  fprintf (stderr, "write output to %s\n", fname);
}
