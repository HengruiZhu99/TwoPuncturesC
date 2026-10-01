/* Native regression: three genuinely independent sequential solves, physical
   sampling, exact Brill-Lindquist data, and an independent Bowen-York tensor.
   The equation residual is recomputed outside the diagnostics interface. */
#include "TwoPunctures.h"

#define REQUIRE(condition) do { if (!(condition)) { \
  fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#condition); \
  exit(EXIT_FAILURE); } } while (0)

static void close_value(double actual, double expected, double tol) {
  REQUIRE(isfinite(actual));
  if (fabs(actual-expected) > tol*(1.0+fabs(expected))) {
    fprintf(stderr,"actual %.17g expected %.17g tolerance %.3g\n",
        actual,expected,tol);
    exit(EXIT_FAILURE);
  }
}

static void configure(double b, double mp, double mm, int na, int nb, int np) {
  TwoPunctures_params_set_default();
  TwoPunctures_params_set_Real("par_b",b);
  TwoPunctures_params_set_Real("par_m_plus",mp);
  TwoPunctures_params_set_Real("par_m_minus",mm);
  TwoPunctures_params_set_Int("npoints_A",na);
  TwoPunctures_params_set_Int("npoints_B",nb);
  TwoPunctures_params_set_Int("npoints_phi",np);
  TwoPunctures_params_set_Int("grid_setup_method",evaluation);
  TwoPunctures_params_set_Int("Newton_maxit",15);
  TwoPunctures_params_set_Real("Newton_tol",1e-11);
}

static double independently_recomputed_residual(ini_data *data) {
  double *residual=malloc(sizeof(double)*data->ntotal);
  derivs *u;
  REQUIRE(residual != NULL);
  allocate_derivs(&u,data->ntotal);
  F_of_v(1,params_get_int("npoints_A"),params_get_int("npoints_B"),
      params_get_int("npoints_phi"),data->v,residual,u);
  double maximum=0.0;
  for (int i=0;i<data->ntotal;++i) {
    REQUIRE(isfinite(residual[i]));
    maximum=fmax(maximum,fabs(residual[i]));
  }
  free_derivs(u);
  free(residual);
  return maximum;
}

/* Abar_ij = 3/(2r^2)[P_i n_j+P_j n_i-(delta_ij-n_i n_j)P.n]
   +3/r^3[(S cross n)_i n_j+(S cross n)_j n_i].
   This function does not call any production Bowen-York routines. */
static void analytic_one_by(const double xyz[3], double location,
    const double P[3], const double S[3], double A[3][3]) {
  const double displacement[3]={xyz[0]-location,xyz[1],xyz[2]};
  const double r=sqrt(displacement[0]*displacement[0]+
      displacement[1]*displacement[1]+displacement[2]*displacement[2]);
  double n[3],dot=0.0;
  for (int i=0;i<3;++i) { n[i]=displacement[i]/r; dot+=P[i]*n[i]; }
  const double cross[3]={S[1]*n[2]-S[2]*n[1],
    S[2]*n[0]-S[0]*n[2],S[0]*n[1]-S[1]*n[0]};
  for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
    A[i][j]+=1.5/(r*r)*(P[i]*n[j]+P[j]*n[i]
        -((i==j ? 1.0 : 0.0)-n[i]*n[j])*dot)
      +3.0/(r*r*r)*(cross[i]*n[j]+cross[j]*n[i]);
  }
}

static void test_rest(double b, double mp, double mm, int na, int nb, int np) {
  configure(b,mp,mm,na,nb,np);
  ini_data *data=TwoPunctures_make_initial_data();
  REQUIRE(data != NULL);
  REQUIRE(data->ntotal == na*nb*np);
  REQUIRE(TwoPunctures_make_initial_data() == NULL); /* one live context */
  TwoPunctures_finalise(NULL); /* does not invalidate the live context */
  const double xyz[15]={1.1,0.2,-0.3, -1.7,0.6,0.1, 0.07,0.03,0.12,
    2.0,-1.5,0.8, 0.4,0.8,-0.3};
  double lapse[5],psi[5],gamma[30],K[30];
  REQUIRE(TwoPunctures_sample_points(data,5,xyz,lapse,psi,gamma,K)==0);
  for (int i=0;i<5;++i) {
    const double x=xyz[3*i],y=xyz[3*i+1],z=xyz[3*i+2];
    const double p=1.0+mp/(2*sqrt((x-b)*(x-b)+y*y+z*z))
      +mm/(2*sqrt((x+b)*(x+b)+y*y+z*z));
    close_value(psi[i],p,3e-14);
    close_value(lapse[i],1.0/(p*p),3e-14);
    for (int q=0;q<6;++q) {
      close_value(gamma[6*i+q],q==0||q==3||q==5 ? pow(p,4) : 0.0,3e-14);
      close_value(K[6*i+q],0.0,0.0);
    }
  }
  double reported,E,masses[2];
  REQUIRE(TwoPunctures_diagnostics(data,&reported,&E,masses)==0);
  const double recomputed=independently_recomputed_residual(data);
  close_value(reported,recomputed,1e-14);
  REQUIRE(recomputed<1e-13);
  close_value(E,mp+mm,1e-14);
  close_value(masses[0],mp+mp*mm/(4*b),1e-14);
  close_value(masses[1],mm+mp*mm/(4*b),1e-14);
  const double puncture[3]={b,0,0};
  REQUIRE(TwoPunctures_sample_points(data,1,puncture,lapse,psi,gamma,K)==-2);
  REQUIRE(TwoPunctures_sample_points(data,0,NULL,NULL,NULL,NULL,NULL)==0);
  REQUIRE(TwoPunctures_sample_points(data,-1,xyz,lapse,psi,gamma,K)==-1);
  printf("BL b=%.3g grid=%dx%dx%d external_Fmax=%.12g ADM=%.12g\n",
      b,na,nb,np,recomputed,E);
  TwoPunctures_finalise(data);
}

static void test_boost_and_spin(void) {
  const double b=0.7;
  configure(b,0.6,0.4,12,10,8);
  const double P[2][3]={{0.02,0.05,-0.015},{-0.02,-0.05,0.015}};
  const double S[2][3]={{0.001,-0.002,0.003},{-0.004,0.001,0.002}};
  char name[40];
  for (int h=0;h<2;++h) for (int a=0;a<3;++a) {
    snprintf(name,sizeof(name),"par_P_%s%d",h==0?"plus":"minus",a+1);
    TwoPunctures_params_set_Real(name,P[h][a]);
    snprintf(name,sizeof(name),"par_S_%s%d",h==0?"plus":"minus",a+1);
    TwoPunctures_params_set_Real(name,S[h][a]);
  }
  ini_data *data=TwoPunctures_make_initial_data();
  REQUIRE(data != NULL);
  REQUIRE(data->ntotal==960);
  const double xyz[12]={1.8,0.3,-0.2, -1.4,0.8,0.1,
    0.2,-0.6,0.4, 0.3,0.7,-0.9};
  double lapse[4],psi[4],gamma[24],K[24];
  REQUIRE(TwoPunctures_sample_points(data,4,xyz,lapse,psi,gamma,K)==0);
  const int row[6]={0,0,0,1,1,2},col[6]={0,1,2,1,2,2};
  double minimum_correction=HUGE_VAL;
  for (int i=0;i<4;++i) {
    double A[3][3]={{0}};
    analytic_one_by(xyz+3*i,b,P[0],S[0],A);
    analytic_one_by(xyz+3*i,-b,P[1],S[1],A);
    const double x=xyz[3*i],y=xyz[3*i+1],z=xyz[3*i+2];
    const double p=1.0+0.6/(2*sqrt((x-b)*(x-b)+y*y+z*z))
      +0.4/(2*sqrt((x+b)*(x+b)+y*y+z*z));
    minimum_correction=fmin(minimum_correction,psi[i]-p);
    close_value(lapse[i],pow(p,-2),3e-14);
    for (int q=0;q<6;++q) {
      close_value(K[6*i+q],A[row[q]][col[q]]/(psi[i]*psi[i]),2e-14);
      close_value(gamma[6*i+q],q==0||q==3||q==5 ? pow(psi[i],4) : 0.0,2e-14);
    }
    close_value(K[6*i]+K[6*i+3]+K[6*i+5],0.0,2e-14);
  }
  REQUIRE(minimum_correction>1e-6); /* includes u, unlike static psi */
  double residual,E,masses[2];
  REQUIRE(TwoPunctures_diagnostics(data,&residual,&E,masses)==0);
  const double external=independently_recomputed_residual(data);
  close_value(residual,external,1e-14);
  REQUIRE(external<1e-9);
  REQUIRE(E>1.0);
  printf("generic BY grid=12x10x8 external_Fmax=%.12g ADM=%.12g min_u=%.12g\n",
      external,E,minimum_correction);

  /* Sampling must agree with the existing Cartesian storage API, including
     its static-factor decomposition and its x/z tensor interchange. Only
     interpolation controls change here; solved physical parameters stay fixed. */
  for (int state=0;state<4;++state) for (int swap=0;swap<2;++swap) {
    params_set_int("conformal_state",state);
    params_set_int("swap_xz",swap);
    params_set_real("center_offset1",0.1);
    params_set_real("center_offset2",-0.05);
    params_set_real("center_offset3",0.03);
    double derivatives[10]={1},oldg[6],oldK[6],old_lapse;
    double x=xyz[0],y=xyz[1],z=xyz[2];
    int zero[3]={0,0,0},one[3]={1,1,1};
    TwoPunctures_Cartesian_interpolation(data,zero,one,one,&x,&y,&z,
        &old_lapse,derivatives,derivatives+1,derivatives+2,derivatives+3,
        derivatives+4,derivatives+5,derivatives+6,derivatives+7,
        derivatives+8,derivatives+9,oldg,oldg+1,oldg+2,oldg+3,oldg+4,oldg+5,
        oldK,oldK+1,oldK+2,oldK+3,oldK+4,oldK+5);
    REQUIRE(TwoPunctures_sample_points(data,1,xyz,lapse,psi,gamma,K)==0);
    for (int q=0;q<6;++q) {
      close_value(gamma[q],oldg[q]*pow(derivatives[0],4),3e-14);
      close_value(K[q],oldK[q],3e-14);
    }
    close_value(lapse[0],old_lapse,3e-14);
  }
  TwoPunctures_finalise(data);
}

int main(void) {
  test_rest(0.3,0.5,0.5,10,10,8);
  test_boost_and_spin();
  test_rest(0.15,0.7,0.2,8,8,6);
  puts("physical API and sequential lifecycle regression passed");
  return EXIT_SUCCESS;
}
