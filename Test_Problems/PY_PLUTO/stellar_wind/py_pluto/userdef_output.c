#include "pluto.h"

double heatcool2();
double sqsqxi;

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k, nv;  



double ***dvdr,Vc[NVAR],g[3];
double ***M,t,v_th,M_max,***T_out;



double rho,p,T,xi,ne,nH,n,mu,lx,tx,r,***gr_out;
double *x1   = grid->x[IDIR],  *x2  = grid->x[JDIR],  *x3  = grid->x[KDIR];


dvdr = GetUserVar("dvdr");
M = GetUserVar("M");
gr_out= GetUserVar("g_r");
T_out = GetUserVar("T");

M_max=g_inputParam[M_rad];


    /*
comp_h_pre  = GetUserVar("ch_pre");
comp_c_pre  = GetUserVar("cc_pre");
xray_h_pre  = GetUserVar("xh_pre");
line_c_pre  = GetUserVar("lc_pre");
brem_c_pre  = GetUserVar("bc_pre");
ne_out	= GetUserVar("ne");
nH_out	= GetUserVar("nh");
T_out     = GetUserVar("T");
xi_out     = GetUserVar("XI");

#if (BODY_FORCE & VECTOR)
gr_out     = GetUserVar("gr");
gt_out     = GetUserVar("gt");
gp_out     = GetUserVar("gp");
#endif


mu=g_inputParam[MU];
lx=g_inputParam[L_x];  //Xray luminosiy
tx=g_inputParam[T_x];  //Xray tenperature 

 */ 
v_th=sqrt(3.)*g_isoSoundSpeed*UNIT_VELOCITY;
  
  
  DOM_LOOP(k,j,i){
  	  r=grid->x[IDIR][i]*UNIT_LENGTH;  //The radius - in real units
	  
      t=0.4*d->Vc[RHO][k][j][i]*UNIT_DENSITY*v_th/fabs(dvdr_array[i]);   
      M[k][j][i]=1./30.*pow(t,-0.7);
      if (M[k][j][i]>M_max) M[k][j][i]=M_max;
      
          
      dvdr[k][j][i]=dvdr_array[i];

	NVAR_LOOP(nv)  Vc[nv]=d->Vc[nv][k][j][i];	
	
	BodyForceVector(Vc, g, x1[i], x2[j], x3[k],i,j,k);
    
#if EOS==ISOTHERMAL
	T_out[k][j][i]=pow(g_isoSoundSpeed*UNIT_VELOCITY,2)*0.6/CONST_Rgas;
#endif    
/* The acceleration vectors are spherical polar */	
	
	gr_out[k][j][i]=g[0];
//	gt_out[k][j][i]=g[1];
//	gp_out[k][j][i]=g[2];
//	
	
//	#endif
	
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





