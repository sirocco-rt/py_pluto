/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Implementation of the isothermal EOS.
                    
  \author A. Mignone (mignone@ph.unito.it)
  \date   March 05, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SoundSpeed2 (const State *p, int beg, int end, int pos, Grid *grid)
/*!
 * Define the square of the sound speed.
 * 
 * \param [in]   p    pointer to a state structure
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 * \param [in]  pos   an integer specifying the spatial position 
 *                    inside the cell (only for spatially-dependent EOS)
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int  i, j,k;  /* -- used as multidimensional indices -- */
  double *x1, *x2, *x3;
  double teff, b, disk_mdot, cent_mass,r_0;
  
  #if PHYSICS == HD || PHYSICS == MHD
   x1 = grid->x[IDIR];
   x2 = grid->x[JDIR];
   x3 = grid->x[KDIR];

   i = g_i; j = g_j; k = g_k;

   disk_mdot = g_inputParam[DISK_MDOT];
   //disk_mdot = 1.e18;
   //cent_mass = 0.6; 
   cent_mass = g_inputParam[CENT_MASS]; 
   r_0=g_inputParam[R_0]/UNIT_LENGTH;
   
   teff = pow(3.*CONST_G*cent_mass*disk_mdot/(8.*CONST_PI*CONST_sigma),0.25);
   teff *= pow(UNIT_LENGTH,-0.75);
 
   if (g_dir == IDIR) {
     double R;
     x1 = (pos == FACE_CENTER ? grid->xr[IDIR] : grid->x[IDIR]);
     //for (i = beg; i <= end; i++) p->a2[i] = g_isoSoundSpeed*g_isoSoundSpeed;
     for (i = beg; i <= end; i++) {
         R = x1[i]*sin(x2[j]);
         if (R<r_0) R = 1.0001*r_0;
         //b = 1.- pow(r_0/R,0.5);
         b = 1.;
         p->a2[i] = CONST_Rgas*teff*pow(R*R*R/b,-0.25)/0.6/(UNIT_VELOCITY*UNIT_VELOCITY);
         //print("i=%i teff %e\n",i, teff*pow(R*R*R/b,-0.25));
     }

   }else if (g_dir == JDIR){
     double R;
     x2 = (pos == FACE_CENTER ? grid->xr[JDIR] : grid->x[JDIR]);
     //for (j = beg; j <= end; j++) p->a2[j] = g_isoSoundSpeed*g_isoSoundSpeed;
     for (j = beg; j <= end; j++) {
         R = x1[i]*sin(x2[j]);
         if (R<r_0) R = 1.0001*r_0;
         //b = 1.- pow(r_0/R,0.5);
         b = 1.;
         p->a2[j] = CONST_Rgas*teff*pow(R*R*R/b,-0.25)/0.6/(UNIT_VELOCITY*UNIT_VELOCITY);
         //print("j=%i teff %e\n",j, teff*R*R*R/b,-0.25));
     }

   }else if (g_dir == KDIR){
     double R;
     x3 = (pos == FACE_CENTER ? grid->xr[KDIR] : grid->x[KDIR]);
     //for (k = beg; k <= end; k++) p->a2[k] = g_isoSoundSpeed*g_isoSoundSpeed;
     for (k = beg; k <= end; k++) {
         R = x1[i]*sin(x2[j]);
         if (R<r_0) R = 1.0001*r_0;
         //b = 1.- pow(r_0/R,0.5);
         b = 1.;
         p->a2[k] = CONST_Rgas*teff*pow(R*R*R/b,-0.25)/0.6/(UNIT_VELOCITY*UNIT_VELOCITY);
         //print("k=%i teff %e\n",k, teff*pow(R*R*R/b,-0.25));
     }

  }
  #else
   printLog ("! SoundSpeed2: not defined for this EoS\n");
   QUIT_PLUTO(1);
  #endif
}

/* ********************************************************************* */
void Enthalpy (double **v, double *h, int beg, int end)
/*!
 * Compute the enthalpy.
 *
 * \param [in]    v   1D array of primitive quantities
 * \param [in]    h   1D array of enthalpy values
 * \param [in]  beg   initial index of computation 
 * \param [in]  end   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  printLog ("! Enthalpy: enthalpy not defined in isothermal EOS\n");
  QUIT_PLUTO(1);
}

/* ********************************************************************* */
void Entropy (double **v, double *s, int beg, int end)
/*!
 * Compute the entropy.
 * 
 * \param [in]    v   1D array of primitive quantities
 * \param [in]    s   1D array of entropy values
 * \param [in]   is   initial index of computation 
 * \param [in]   ie   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  printLog ("Entropy: entropy not defined for isothermal EOS\n");
  QUIT_PLUTO(1);
}
