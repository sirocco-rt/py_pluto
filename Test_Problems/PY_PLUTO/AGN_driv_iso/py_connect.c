/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief  Read in a py_heatcool file and compare with known heating and cooling rates

n constant
  
 
  \authors Nick H
  \date    Dec 3 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define LINELENGTH 800



/* ********************************************************************* */
void read_py_heatcool (Data *d, Grid *grid,int flag)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
	int i,j,k;
	FILE *fopen (), *fptr;
	char *fgets (), aline[LINELENGTH];
	int ii,jj,nwords;
	double dens,comp_h_pre,comp_c_pre,xray_h_pre,brem_c_pre,line_c_pre;
	double rcen,thetacen;
	int icount;
	#if (BODY_FORCE & VECTOR)
	double gr,gt,gp;
	#endif
	

	if (flag==0) //Initialising first time round
    {
    	printf ("Entering py_heatcool - first time round\n");
        
    DOM_LOOP(k,j,i){
        
#if COOLING == BLONDIN
		d->comp_h_pre[k][j][i]=1.0;
		d->comp_c_pre[k][j][i]=1.0;
		d->xray_h_pre[k][j][i]=1.0;
		d->line_c_pre[k][j][i]=1.0;
		d->brem_c_pre[k][j][i]=1.0;
#endif	
	}
	printf ("Initialised\n");
}
	else //We will be reading in some updated numbers
	{
		fptr = fopen ("prefactors.dat", "r");
		if (fptr==NULL)
		{
		    DOM_LOOP(k,j,i){
#if COOLING == BLONDIN                
				d->comp_h_pre[k][j][i]=1.0;
				d->comp_c_pre[k][j][i]=1.0;
				d->xray_h_pre[k][j][i]=1.0;
				d->line_c_pre[k][j][i]=1.0;
				d->brem_c_pre[k][j][i]=1.0;
#endif                		
			}
			printf ("NO prefactor file\n");
		}
		else
		{
	    DOM_LOOP(k,j,i){ //Initialise
#if COOLING == BLONDIN                            
			d->comp_h_pre[k][j][i]=-1.0;
			d->comp_c_pre[k][j][i]=-1.0;
			d->xray_h_pre[k][j][i]=-1.0;
			d->line_c_pre[k][j][i]=-1.0;
			d->brem_c_pre[k][j][i]=-1.0;
#endif            
		    }
        if (fgets (aline, LINELENGTH, fptr) != NULL) {
            
        } else {
            printf("Error reading aline\n");
            exit(0);
        }
		icount=0;
		while (fgets (aline, LINELENGTH, fptr) != NULL)	
		{		
			if 	((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le %le %le %le %le %le", &ii, &rcen, &jj, &thetacen, &dens,
				&comp_h_pre, &comp_c_pre, &xray_h_pre, &brem_c_pre, &line_c_pre, &gr, &gt, &gp)) == 13)
				{
					DOM_LOOP(k,j,i)
					{
						if (fabs(((rcen/UNIT_LENGTH)-grid->x[IDIR][i])/(rcen/UNIT_LENGTH))<1e-5 && fabs((thetacen-grid->x[JDIR][j])/thetacen)<1e-5)
						{
							if ((d->Vc[RHO][k][j][i]*UNIT_DENSITY/dens)-1.>1e-6)
							{
								printf ("Density mismatch in cell i=%i j=%i old=%e new=%e\n",i,j,d->Vc[RHO][k][j][i]*UNIT_DENSITY,dens);
							}
							icount++;
							#if COOLING == BLONDIN                                                    
							d->comp_h_pre[k][j][i]=comp_h_pre;
							d->comp_c_pre[k][j][i]=comp_c_pre;
							d->xray_h_pre[k][j][i]=xray_h_pre;
							d->line_c_pre[k][j][i]=line_c_pre;
							d->brem_c_pre[k][j][i]=brem_c_pre;
							#endif
							}
						}
					}
				else
				{
					printf ("Prefactor file incorrectly formatted nwords=%i %s\n",nwords,aline);
					exit(0);
				}
			}
#if COOLING == BLONDIN                                                            
	    DOM_LOOP(k,j,i){ //Test
			if (d->comp_h_pre[k][j][i]<0.0)
			{
				printf ("comp_h_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->comp_h_pre[k][j][i]);
			}
			if (d->comp_c_pre[k][j][i]<0.0)
			{
				printf ("comp_c_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->comp_c_pre[k][j][i]);
			}
			if (d->xray_h_pre[k][j][i]<0.0)
			{
				printf ("xray_h_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->xray_h_pre[k][j][i]);
			}
			if (d->line_c_pre[k][j][i]<0.0)
			{
				printf ("line_c_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->line_c_pre[k][j][i]);
			}
			if (d->brem_c_pre[k][j][i]<0.0)
			{
				printf ("brem_c_prefactor<0.0 for i=%i j=%i comp_h_pre=%e\n",i,j,d->brem_c_pre[k][j][i]);
			}
		}
#endif
        
    
	
		printf ("Read in %i prefectors\n",icount);
	}
}
}



/*This routine is used when we are running isothermal runs with driving*/
#if EOS==ISOTHERMAL && PY_CONNECT
/* ********************************************************************* */
void read_py_iso_temp (Data *d, Grid *grid,int flag)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
	int i,j,k;
	FILE *fopen (), *fptr;
	char *fgets (), aline[LINELENGTH];
	int ii,jj,nwords;
	double dens,temp,nh,ne,t_opt,t_UV,t_Xray,r,x1in,x2in;
	int icount;
	#if (BODY_FORCE & VECTOR)
	double gx,gy,gz;
	#endif
    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    
    double r1=2e12/UNIT_LENGTH;
    double t0=50000;
    double t1=35000; 
    double dtdr; 
    double tol=1e-6;
	
    
    dtdr=(t1-t0)/(r1-grid->x[IDIR][0]); 
    
    
    temp=g_inputParam[T_ISO];
	
    if (flag==1)
    {
        printf ("We are restarting so expecting to read in temperatures\n");
        fptr = fopen ("py_pcon_data.dat", "r");
        if (fptr==NULL)
        {
            printf ("NO pcon file\n");
    	}
        if (fgets (aline, LINELENGTH, fptr) != NULL){
            
        } else {
            printf("Error in reading aline\n");
            exit(0);
        }
		
		icount=0;
    	while (fgets (aline, LINELENGTH, fptr) != NULL)	
		{
    	    if ((nwords = sscanf (aline, "%d %d %le %le %le %le %le %le %le %le %le", &ii, &jj, &x1in, &x2in, &temp, &dens, &nh, &ne, &t_opt, &t_UV, &t_Xray)) == 11)
			//Now we have to find what cell it relates to - we look for the geometry
			DOM_LOOP(k,j,i)
			{
				if (fabs(1.-(x1in/UNIT_LENGTH/x1[i]))<tol && fabs(1.-(x2in/x2[j]))<tol)
				{
                    py_temp[k][j][i]=temp;
					icount++;
				}
			}
            else
            {
                printf ("Temerature file incorrectly formatted\n");
            }
		}
		printf ("Read temperatures for  %i cells\n",icount);
	    fclose(fptr);
		
		
        TOT_LOOP(k,j,i) //Fill in the ends of the temperature array
        {
            if (i<IBEG)
            {
                py_temp[k][j][i]=py_temp[k][j][IBEG];
            }
            else if (i>IEND)
            {
                py_temp[k][j][i]=py_temp[k][j][IEND];
            }
        }
    }
    else
    {
    	printf ("First time round so setting temperatures to fixed T\n");
        TOT_LOOP(k,j,i)
        {
            if (grid->x[IDIR][i] < r1)
            {
                py_temp[k][j][i]=t0+(grid->x[IDIR][i]-grid->x[IDIR][0])*dtdr;
            }
            else
            {
                py_temp[k][j][i]=t1;//For the initial run - we have no driving                	
            }
            py_temp[k][j][i]=temp;
        } 
               	
    }
    
  
    	
}

#endif

/*This routine is used when we are running isothermal runs with driving*/

/* ********************************************************************* */
void read_py_rad_driv (Data *d, Grid *grid,int flag)
/*!
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
	int i,j,k;
	FILE *fopen (), *fptr;
	char *fgets (), aline[LINELENGTH];
	int ii,jj,nwords;
	double rcen,thetacen,dens,gx,gy,gz;
	int icount;

	printf ("Reading in accelerations from python\n");
    if (flag==1)
    {
        fptr = fopen ("py_accelerations.dat", "r");
        if (fptr==NULL)
        {
            printf ("NO accelerations file\n");
    	}
        if (fgets (aline, LINELENGTH, fptr) != NULL){
            
        } else {
            printf("Error reading aline\n");
            exit(0);
        }
    	icount=0;
    	while (fgets (aline, LINELENGTH, fptr) != NULL)	
    	{		
    	    if ((nwords = sscanf (aline, "%d %le %d %le %le %le %le %le",  &ii, &rcen, &jj, &thetacen, &dens, &gx, &gy, &gz)) == 8)
    	    {
//                if (fabs((dens-d->Vc[RHO][k][jj+JBEG][ii+IBEG]*UNIT_DENSITY)/dens) > 1e-6)
//                    printf ("NOT SAME  %e %e %e %e\n",thetacen,dens,gx,d->Vc[RHO][k][jj+JBEG][ii+IBEG]*UNIT_DENSITY);
//                else
//                {
                    icount++;
    				g_rad[0][k][jj+JBEG][ii+IBEG]=gx/UNIT_ACCELERATION;
    				g_rad[1][k][jj+JBEG][ii+IBEG]=gy/UNIT_ACCELERATION;
    				g_rad[2][k][jj+JBEG][ii+IBEG]=gz/UNIT_ACCELERATION;
//                }
                
            }
            else
            {
                printf ("Incorrectly formatted\n");
            }
    	}
        printf ("matched %i accelerations\n",icount);
    }
    else
    {
        DOM_LOOP(k,j,i){
        	printf ("WE ARE HERE - first time round\n");
    		g_rad[0][k][j][i]=12345678.; //For the initial run - we have no driving - set a flag
    		g_rad[1][k][j][i]=12345678.;
    		g_rad[2][k][j][i]=12345678.;		
    	}        
    }	
}



void read_py_fluxes (Data *d, Grid *grid)
/*
 * \param [in,out]  VV    a pointer to the PLUTO 3D data array containing
 *                        pimitive variables.
 * \param [in]      dt    the current integration time step
 * \param [in]      Dts   a pointer to the Time_Step structure
 * \param [in]      grid  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
	int i,j;
	long int k;
	FILE *fopen (), *fptr;
	char *fgets (), aline[LINELENGTH];
	long int ii,jj,kk,nwords;
	double x1cen,x2cen,xin,zin,x1in,x2in;
    double kradin,alpharadin,rhoin;
    double krad,alpharad;
    double opt[4],UV[4],Xray[4];
	double temp;
	int iflux;
	int icount,iaxis,match;
	double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];
    double *x1l = grid->xl[IDIR];
    double *x2l = grid->xl[JDIR];
    double *x3l = grid->xl[KDIR];
	double *x1r = grid->xr[IDIR];
    double *x2r = grid->xr[JDIR];
    double *x3r = grid->xr[KDIR];
	double rcen,thetacen,dens,gx,gy,gz;
    double tol;
    kk=0;
    krad=g_inputParam[KRAD];               
    alpharad=g_inputParam[ALPHARAD];
	
	
	
	printf ("Arrived in read_py_fluxes IBEG=%li JBEG=%li\n",IBEG,JBEG);
	printf ("Reading in fluxes from python\n");
	tol = 1.e-2;
	for (iaxis=0;iaxis<3;iaxis++) //We loop over three axes - x,y,z in seperate files
	{
		if (iaxis==0) fptr = fopen ("directional_flux_theta.dat", "r");
		else if (iaxis==1) fptr = fopen ("directional_flux_phi.dat", "r");
		else if (iaxis==2) fptr = fopen ("directional_flux_r.dat", "r");
	    if (fptr==NULL)
	    {
	        printf ("No flux file\n");
		}
        if (fgets (aline, LINELENGTH, fptr) != NULL){
            fgets(aline, LINELENGTH, fptr);
        } else {
            printf("Error reading aline\n");
            exit(0);
        }
		nwords = sscanf (aline, "%*s %*s %ld ",  &ii);
		if (nwords == 1)
		{
			if (iaxis==0)
			{			
				NFLUX_ANGLES=ii;
				printf ("We have %i angular bins for flux\n",NFLUX_ANGLES);	
				printf ("Reading in X fluxes\n");
		    	flux_x_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);  //This is a global variable because rad force doesnt have access to data	
				printf ("declared an array of size %li %li %li\n",NX3_TOT, NX2_TOT, NX1_TOT);
			}
			else if (iaxis==1)
			{
				if (ii!=NFLUX_ANGLES)
				{
					printf ("Y-flux file doesnt agree in NFLUX_ANGLES CRASH!!\n");
					exit(0);
				}
				else
				{	
			    	flux_y_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);  //This is a global variable because rad force doesnt have access to data
					printf ("Reading in Y fluxes\n");

				}			
			}
			else if (iaxis==2)
			{
				if (ii!=NFLUX_ANGLES)
				{
					printf ("Z-flux file doesnt agree in NFLUX_ANGLES CRASH!!\n");
					exit(0);
				}
				else
				{
			    	flux_z_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);  //This is a global variable because rad force doesnt have access to data
					printf ("Reading in Z fluxes\n");
					
				}			
			}
		}
		else
		{
			printf ("Flux header improperly formatted\n");
		}
		icount=0;
		match=0;
        //printf("hello\n");
		while (fscanf(fptr,"%ld ",&ii) !=EOF) //Read the first element of a line, and check for EOF
		{
            if (fscanf(fptr,"%ld %*ld %le %le",&jj,&xin, &zin) == 3) //Get the rest of the geometry stuff for this line
            {
                //Now we have to find what cell it relates to - we look for the geometry
                DOM_LOOP(k,j,i)
                {
					x1cen = 0.5*(x1l[i]+x1r[i])*sin(0.5*(x2l[j]+x2r[j]));
					x2cen = 0.5*(x1l[i]+x1r[i])*cos(0.5*(x2l[j]+x2r[j]));

					//printf("x2[i+1]=%e, x2[i]=%e\n",x2[i+1],x2[i]);

                    if (fabs(1.-(xin/UNIT_LENGTH/x1cen))<tol && fabs(1.-(zin/UNIT_LENGTH/x2cen))<tol)
                    {
                        //printf("hello\n");
                        //printf ("found a match %e %e %i %e %e %i\n",xin/UNIT_LENGTH,x1cen,i,zin/UNIT_LENGTH,x2cen,j);
                        for (iflux=0;iflux<NFLUX_ANGLES;iflux++)
                        {
                            if (fscanf (fptr,"%le",&temp) == 1) //read the flux in one angular bin at a time
                            {
                                if (iaxis==0)
                                    flux_x_UV[iflux][k][j][i]=temp;
                                    //printf("flux_x_UV=%e\n",temp);
                                if (iaxis==1)
                                    flux_y_UV[iflux][k][j][i]=temp;
                                if (iaxis==2)
                                    flux_z_UV[iflux][k][j][i]=temp;
                            }
                            else
                            {
                                printf("Error in reading flux file\n");
                                exit(0);
                            }
                        }
                        match=1;
                        icount++;
                    } //else {
                        //printf("x1in=%e, x2in=%e, fabs_r=%e,fabs_theta=%e\n",x1in, x2in, fabs(1.-(x1in/UNIT_LENGTH/x1[i])),fabs(1.-(x2in/x2[j])));
                    //}
                }
                if (match==0)
                {
                    for (iflux=0;iflux<NFLUX_ANGLES;iflux++)
                    {
                        if (fscanf (fptr,"%le",&temp) ==1)
                        {  //read the flux in one angular bin at a time
                        }
                        else
                        {
                            printf("Error in reading flux file because no match\n");
                            exit(0);
                        }
                    }
                }
                match=0;
            } else
            {
                printf("Error in reading flux file\n");
                exit(0);
            }
        }
		printf ("Read %i fluxes for %i cells\n",NFLUX_ANGLES,icount);
	    fclose(fptr);
	}	
	
	flux_r_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
	flux_t_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
	flux_p_UV = ARRAY_4D(NFLUX_ANGLES,NX3_TOT, NX2_TOT, NX1_TOT, double);
	
//Make fluxes in r,theta,phi directions as well	
	
	for (iflux=0;iflux<NFLUX_ANGLES;iflux++)
	{ 
		DOM_LOOP(k,j,i)
		{
			//flux_r_UV[iflux][k][j][i]=flux_x_UV[iflux][k][j][i]*sin(x2[j])+flux_z_UV[iflux][k][j][i]*cos(x2[j]);
			flux_r_UV[iflux][k][j][i]=flux_z_UV[iflux][k][j][i];
			//flux_t_UV[iflux][k][j][i]=flux_x_UV[iflux][k][j][i]*cos(x2[j])-flux_z_UV[iflux][k][j][i]*sin(x2[j]);
			flux_t_UV[iflux][k][j][i]=flux_x_UV[iflux][k][j][i];
			flux_p_UV[iflux][k][j][i]=flux_y_UV[iflux][k][j][i];		
		}
	}
	
	
	k=0; //Need to reset this number to zero after the dom loop!
	
	tol = 1.e-6;
    if (krad==999 && alpharad==999)
    {
		printf ("Importing a force multiplier fit file\n");
		fptr=fopen ("M_UV_data.dat", "r");
	    if (fptr==NULL)
	    {
	        printf ("No force multiplier file\n");
		}
        if (fgets (aline, LINELENGTH, fptr) != NULL){
            
        } else {
            printf("Error reading aline\n");
            exit(0);
        }
		nwords = sscanf (aline, "%*s %ld ",  &ii);
		MPOINTS=ii;
		printf ("There are %i points in the force multiplier fits\n",MPOINTS);
		M_UV_fit=ARRAY_4D(MPOINTS,NX3_TOT, NX2_TOT, NX1_TOT, double);
		
        if (fscanf(fptr,"%*s ")==0){
                //Get the rest of the geometry stuff for this line
        } else {
            printf("Error in reading geometry from force multiplier file\n");
            exit(0);
        }
		t_fit=calloc(MPOINTS, sizeof(double));
		
        for (iflux=0;iflux<MPOINTS;iflux++)
		{
			if (fscanf (fptr,"%le",&temp) == 1) //read the values of t for which the force multiplier is tabulates
            {
                t_fit[iflux]=log10(temp);
                //printf ("BOOM %e\n",t_fit[iflux]);
            } else
            {
                printf("Error in reading force multiplier file\n");
                exit(0);
            }
		}
		icount=0;
		match=0;
		while (fscanf(fptr,"%ld ",&ii) !=EOF) //Read the first element of a line, and check for EOF
		{
			if (fscanf(fptr,"%ld %le %le",&jj,&x1in,  &x2in) == 3) //Get the rest of the geometry stuff for this line
            {
                //			printf ("Processing cell %i %i\n",ii,jj);
                DOM_LOOP(k,j,i)
                {
                    if (fabs(1.-(x1in/UNIT_LENGTH/x1[i]))<tol && fabs(1.-(x2in/x2[j]))<tol)
                    {
                        //printf ("found a match %e %e %i %e %e %i\n",x1in,x1[i]*UNIT_LENGTH,i,x2in,x2[j],j);
                        for (iflux=0;iflux<MPOINTS;iflux++)
                        {
                            if (fscanf (fptr,"%le",&temp) == 1){ //read the flux in one angular bin at a time
                                M_UV_fit[iflux][k][j][i]=log10(temp);
                            } else {
                                printf("Error in reading force multiplier file\n");
                                exit(0);
                            }
                        }
                        match=1;
                        icount++;
                    }
                }
                if (match==0)
                {
                    for (iflux=0;iflux<MPOINTS;iflux++)
                    {
                        if (fscanf (fptr,"%le",&temp) == 1){ //read the flux in one angular bin at a time
                        } else {
                            printf("Error in reading force multiplier file\n");
                            exit(0);
                        }
                    }
                }
                match=0;
            } else
            {
                printf("Error in reading force multiplier file\n");
                exit(0);
            }
            
		}	    
		fclose(fptr);
		printf ("Read %i points to M vs t fits for %i cells\n",MPOINTS,icount);
		
    }
    else
    {
        printf ("Using k=%e alpha=%e for all cells\n",krad,alpharad);
    }
        printf ("Not using accelerations from python\n");

DOM_LOOP(k,j,i)
                {
                        g_rad[0][k][j][i]=0.0;
                        g_rad[1][k][j][i]=0.0;
                        g_rad[2][k][j][i]=0.0;
                } 

}




