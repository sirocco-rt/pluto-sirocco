#include "cak_v2.h"


#define f_uv_lower 7.5e14
#define f_uv_upper 3e16

#define NPOINTS 32 //This is the number of t (and hence M) bins
#define tmin -15
#define tmax 5

int my_rank;

double scale=1.0;

double bb_temp;
double bb_w;

/* Main - top level routine. */
int main(int argc, char** argv)
{

  FILE *output;
  
  double line_nu;
  double A_ul, B_ul, B_lu;
  double kappa_l, delta_doppler, flux_factor, test;
  int upper, lower;

  int icell,iband,itrans,iion,ilev;
  int n_used_lines;
  int v_uv_x; //An index used to decide if we are incremeneting the vis(0), uv(1)or xray(2) band 

//  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel
  int my_nmin,my_nmax;
  int np_mpi_global,num_mpi_cells,num_mpi_extra,rank_global;
  int i,ii,it;
  double *part;
  double t_array[NPOINTS];
  
  double *M_array_transmit;
  double *M_array_transmit2;
  
  double dt;


#ifdef MPI_ON
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif
np_mpi_global = np_mpi;
rank_global=my_rank;
  g_usedMemory=0.0;
  


  for (i=0;i<NPOINTS;i++)
  {
	  dt=(tmax-tmin)/(NPOINTS-1.0);
	  t_array[i]=pow(10,tmin+i*dt);
	  if (my_rank==0)
	  {
	  printf ("BOOM %e %e\n",t_array[i],dt);
  	  }
  }

  read_ionfs();  //Read in the ion fractions and also the     
    
  read_cont();  //Read in the continuum
  
  read_line_data(); //Read in all the line data

	M_UV_array=ARRAY_2D(ncells,NPOINTS,double); //We compute M for NPOINTS values of t
	part=calloc(nions,sizeof(double));
	  
	if (my_rank==0) printf ("Set up working arrays - currently using %f Mb\n",g_usedMemory/1e6);
	if (my_rank==0) printf ("Calculating\n");
	
    
    my_nmin = 0;
    my_nmax = ncells;
#ifdef MPI_ON
  num_mpi_cells = floor (ncells / np_mpi_global);
  num_mpi_extra = ncells - (np_mpi_global * num_mpi_cells);
  
  if (my_rank==0)  printf ("BLAH num_mpi_cells=%i num_mpi_extra=%i\n",num_mpi_cells,num_mpi_extra);

  /* this section distributes the remainder over the threads if the cells
     do not divide evenly by thread */
  if (rank_global < num_mpi_extra)
  {
    my_nmin = rank_global * (num_mpi_cells + 1);
    my_nmax = (rank_global + 1) * (num_mpi_cells + 1);
  }
  else
  {
    my_nmin = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra) * (num_mpi_cells);
    my_nmax = num_mpi_extra * (num_mpi_cells + 1) + (rank_global - num_mpi_extra + 1) * (num_mpi_cells);
  }
#endif
    
    
    
    
    
  if (my_rank==0) printf ("Thread 0 doing cells %i to %i\n",my_nmin,my_nmax);
    
//	for (icell=0; icell<ncells; icell++)  //We are gping to loop over all cells
    for (icell = my_nmin; icell < my_nmax; icell++)
	{
		n_used_lines=0;
/* first thing we need to do is set up the level populations for this cell*/
		for (iion=0;iion<nions;iion++)  //loop over all ions
		{            
			part[iion]=0.0;
			for (ilev=0;ilev<ion_info_nlevels[iion];ilev++) //loop over all the levels for this ion
			{
//                if (ilev>NLEVELS) printf ("BOOM\n");
				lev_pops[iion][ilev]=lev_weight[iion][ilev]*exp1(-1. * lev_energy[iion][ilev] *EV/KB/T_e[icell]); //Level population - based on boltzmann distribution
				part[iion] += lev_pops[iion][ilev]; 					
			}
			part[iion] = 1./part[iion];	
			for (ilev=0;ilev<ion_info_nlevels[iion];ilev++) //loop over all the levels for this ion
			{
				lev_pops[iion][ilev]*=part[iion]*ion_fracs[icell][iion]; //Level population - based on boltzmann distribution
			}			
		}					
/*We now have level populations and partition function for all the ions in this cell*/ 		
		for (itrans=0;itrans<nline_tot;itrans++) //Loop over all lines
		{
			if (trans_freq[itrans] > 0.0) 
			{
				iion = trans_ion[itrans]; //This gets the link into the various ion related arrays
//                if (iion>NIONS) printf ("BOOM\n");
				if (ion_fracs[icell][iion]>0.0)
				{
					if (trans_freq[itrans]<f_uv_lower) v_uv_x=0;
					else if (trans_freq[itrans]>f_uv_upper) v_uv_x=2;
					else v_uv_x=1;
					n_used_lines++;
					A_ul = trans_A_ul[itrans];
					B_ul = trans_B_ul[itrans];
					B_lu = trans_B_lu[itrans];
					lower= trans_lower[itrans];
					upper= trans_upper[itrans];
					line_nu= trans_freq[itrans];
					kappa_l = (B_lu*lev_pops[iion][lower] - B_ul*lev_pops[iion][upper]) * param1[icell]; 
					delta_doppler = line_nu * v_thermal[icell]/C;
					if (J[icell]==0.0 || trans_iband[itrans]<0)
						flux_factor=0.0; //There is no flux, and so no force multiplier
					else
					{                       
                    	if (model_type==IONMODE_MATRIX_SPECTRALMODEL)  // matrix solver spectral model based on power laws
                    	{
						flux_factor = model_jnu(line_nu, trans_lfreq[itrans],icell,trans_iband[itrans],0)/J[icell];	
                        }
                        else if (model_type==IONMODE_MATRIX_BB || model_type == IONMODE_ML93)
                        {
                            flux_factor = bb2(trans_freq[itrans],W[icell],t_r[icell])/J[icell];	         
                        }  
                        else
                        {
                            printf ("UNKNOWN MODEL TYPE - ABORTING\n");
                            exit(0);
                        }
						for (it=0;it<NPOINTS;it++)
	                    { 
	    					if ( (test = kappa_l*t_array[it]/sigma_e[icell]) > 1.e-6) //If we have a large opacity
	    					{	
	    						M_UV_array[icell][it] += delta_doppler * flux_factor * (1. - exp1(-1.*test)) / t_array[it];  //increment the fore multiplier for this t
//	                            by_ion[icell][iion]+=  delta_doppler * flux_factor * (1. - exp1(-1.*test)) / t_array[it];
                            
	    						if (!isfinite(M_UV_array[icell][v_uv_x]))
	    						{
	    							printf("Non-finite M1. Abort.\n");
	    							printf("%i %i %i %g %g %g %g\n", icell,cell_i[icell],cell_j[icell],delta_doppler,flux_factor,exp1(-1.*test), t_array[it]);
	    							exit(0);
	    						}				
	    					}
	    					else
	    					{
	    						M_UV_array[icell][it] += delta_doppler * flux_factor * kappa_l / sigma_e[icell];	
//	                            by_ion[icell][iion]+=  delta_doppler * flux_factor * kappa_l / sigma_e[icell];
                            					
	    						if (!isfinite(M_UV_array[icell][v_uv_x]))
	    						{
	    							printf("Non-finite M2. Abort.\n");
	    							printf("cell=%i ddop=%g ff=%g J=%e kappa_l=%g line_nu=%g\n",icell, delta_doppler,flux_factor,J[icell],kappa_l,line_nu);
	                                model_jnu(line_nu, trans_lfreq[itrans],icell,trans_iband[itrans],1);
	    							exit(0);
	    						}
	    					}
						} //End of loop over t
	            	} //End of check for non zero flux factor 
				} //End of loop to only do computations if ion frac is greater than zero							
			} //End of loop to compute M if the trans freq is greater than zero 
			else
			{					
				printf ("Dodgy line %i ion %i lev %i-%i  freq %e ignoring\n",itrans,trans_ion[itrans],trans_upper[itrans],trans_lower[itrans],trans_freq[itrans]);
			}					
		}	//End of the loop over all lines
	}//End of the loop over all cells

	if (my_rank==0) printf ("Done all the calculations\n");
    
    

    
    
#ifdef MPI_ON                   // these routines should only be called anyway in parallel but we need these to compile
    MPI_Barrier (MPI_COMM_WORLD);
    
	M_array_transmit=calloc(ncells*NPOINTS,sizeof(double));
	M_array_transmit2=calloc(ncells*NPOINTS,sizeof(double));
    
	
	
	
    for (icell = my_nmin; icell < my_nmax; icell++)
    {
		for (it=0;it<NPOINTS;it++)
		{
			M_array_transmit[icell*NPOINTS+it]=M_UV_array[icell][it];
		}
    }
    MPI_Barrier (MPI_COMM_WORLD);    
    MPI_Reduce (M_array_transmit, M_array_transmit2, ncells*NPOINTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (M_array_transmit2, ncells*NPOINTS, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
    MPI_Barrier (MPI_COMM_WORLD);
    for (icell=0;icell<ncells;icell++)
    {
		for (it=0;it<NPOINTS;it++)
		{
        	M_UV_array[icell][it]=M_array_transmit2[icell*NPOINTS+it];
		}
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    free(M_array_transmit);
    free(M_array_transmit2);
	if (my_rank==0) printf ("Communicated all the arrays\n");
        
#endif
    
    
#ifdef MPI_ON
  if (rank_global == 0)
  {
#endif
      
      printf ("outputting the file\n");

    if ((output = fopen("M_UV_data.dat", "w")) == NULL)
      {
	printf("Cannot open M__UV_data.dat\n");
	abort();
      }
	fprintf(output, "MPOINTS %i\n",NPOINTS);
	fprintf(output, "t ");
	for (it=0;it<NPOINTS;it++)
	{
		fprintf (output,"%e ",t_array[it]);
	}
	fprintf (output,"\n ");
	
			  	  
	for (icell = 0; icell < ncells; icell++)
	{
		fprintf(output, "%i %i ", cell_i[icell], cell_j[icell]);
		for (it=0;it<NPOINTS;it++)
		{
			fprintf (output,"%e ",M_UV_array[icell][it]);
		}
		fprintf(output,"\n");
	}
    fclose(output);
    
    printf ("output the file\n");
    
/*    
    if ((output = fopen("M_data_by_ion.dat", "w")) == NULL)
        {
        printf("Cannot open M_data.dat\n");
        abort();
        }
        
    fprintf (output, "nions %i\n", nions);
        for (i = 0; i < nions; i++)
        {
          fprintf (output, "ion %i dummy %i %i\n", i, ion_info_z[i], ion_info_state[i]);
        }
        fprintf (output, "nplasma %i\n", ncells);
        
    for (icell = 0; icell < ncells; icell++)
    	{
    		fprintf(output, "%i %i ", cell_i[icell], cell_j[icell]);
            
                for (ii = 0; ii < nions; ii++)
                {
          fprintf (output, "%e ", by_ion[icell][ii]);
      }
        fprintf (output, "\n");
    }
    
*/
    
    
#ifdef MPI_ON
  }
#endif
    
    
  /*
  printf ("1");
  
//  free(M_array);
  free(part);
  free(ion_info_z);
  free(ion_info_state);
  free(ion_info_nlevels);
  free(ion_info_nlines);
  free(ion_info_xi);	
  free(cell_i); 
  free(cell_j); 
  free(T_e); 
  free(nnh);
  free(nne);
  free(rho);
  free(v_thermal);
  free(sigma_e);
  free(param1);
//  free(t);
  printf ("2");
//  free(ion_fracs);  
//  free(f1);
//  free(f2);
//  free(model);	
//  free(pl_w);
//  free(pl_alpha);
//  free(exp_w);
//  free(exp_temp);
  
  free(mod_fmin);
  free(mod_fmax);
  free(J);
  
//  free(lev_energy);
//  free(lev_weight);
//  free(lev_pops);
  
  free(trans_ion);
  free(trans_lower);
  free(trans_upper);
  free(trans_iband);
  printf ("3");


  free(trans_freq);
  free(trans_lfreq);
  free(trans_A_ul);
  free(trans_B_ul);
  free(trans_B_lu);
  printf ("4");*/

#ifdef MPI_ON
  MPI_Finalize ();
#endif

}


/**************************************************************************************/
//																					
// This bit of code reads in the ion details.	
//
// It populates the array ion_fracs[NELEMS][NIONS] which are the abundances of ions
//
/**************************************************************************************/


int read_ionfs()
{
  int nelem,i,j;
  FILE *abund,*pcon;
  char filename[100];
  char junk[LINELENGTH];
  double dum1,dum2,dum3,dum4,dum5,dum6,dum7;
  int idum1,idum2;
  char *found;
  
/* A loop over the elements - NB at present we need *all* the elements in place - in files which comtain the number of the element as the filename*/  


  if ((abund=fopen("py_ion_data.dat","r")) ==NULL)
  {
	  printf ("Cannot open py_ion_data.dat\n");
	  abort();
  }
  
	fscanf(abund, "%*s %d", &idum1);
	nions=idum1;
	if (my_rank==0) printf ("We have %i ions\n",nions);
	
  	fgets(junk, 1000, abund); 
	
	
	ion_info_z=calloc(nions,sizeof(int));
	ion_info_state=calloc(nions,sizeof(int));
	ion_info_nlevels=calloc(nions,sizeof(int));
	ion_info_nlines=calloc(nions,sizeof(int));
	ion_info_xi=calloc(nions,sizeof(double));

	g_usedMemory+=4.*nions*sizeof(int);
	g_usedMemory+=nions*sizeof(double);

	
	for (i=0;i<nions;i++)
	{
		if (fgets (junk, LINELENGTH, abund) == NULL)
        {
         	printf ("Error reading ion data %i of %i\n",i,nions);
         	abort();
        }
		sscanf (junk, "%*s %*d %*s %d %d", &idum1, &idum2);
		ion_info_z[i]=idum1;
		ion_info_state[i]=idum2;
	}
	
		
	fscanf(abund, "%*s %d", &idum1);
	ncells=idum1;	
	if (my_rank==0) printf ("We have %i cells\n",ncells);
    
    
	ion_fracs=ARRAY_2D(ncells,nions,double); //Set up the ion fracs array
	by_ion=ARRAY_2D(ncells,nions,double); //Set up the ion fracs array
		
	for (i=0;i<ncells;i++)
	{
		fscanf(abund, " %d", &idum1);
		fscanf(abund, " %d", &idum2);
		for (j=0;j<nions;j++)
		{
			fscanf(abund, "%le",&dum1);
			ion_fracs[i][j]=dum1;
		}
	}
    fclose(abund);  //Close the file
	
	if (my_rank==0) printf ("Read ion data - currently using %f Mb\n",g_usedMemory/1e6);
	
	
	
    if ((pcon=fopen("py_pcon_data.dat","r")) ==NULL)
    {
  	  printf ("Cannot open py_pcon_data.dat\n");
  	  abort();
    }	
	
	if (fgets (junk, LINELENGTH, pcon) == NULL)
    {
     	printf ("Error reading first line of py_pcon_data.dat\n");
     	abort();
    }
		
	sscanf(junk, "%*s %d", &idum1);
	
	
	if (idum1!=ncells)
	{
		printf ("py_pcon_data has different number of cells - aborting");
		abort();
	}

	cell_j=calloc(ncells,sizeof(int)); 
	cell_i=calloc(ncells,sizeof(int));

	T_e=calloc(ncells,sizeof(double)); 
	nnh=calloc(ncells,sizeof(double));
	nne=calloc(ncells,sizeof(double));
	rho=calloc(ncells,sizeof(double));
	v_thermal=calloc(ncells,sizeof(double));
	sigma_e=calloc(ncells,sizeof(double));
	param1=calloc(ncells,sizeof(double));
	
	g_usedMemory+=7.*ncells*sizeof(double);
	g_usedMemory+=2.*ncells*sizeof(int);
	
	
	t=ARRAY_2D(ncells,3,double); //Three bands, vis, UV and Xray
	
	
	for (i=0;i<ncells;i++)
	{
		if (fgets (junk, LINELENGTH, pcon) == NULL)
        {
         	printf ("Error reading pcon data %i of %i\n",i,ncells);
         	abort();
        }
		sscanf (junk, "%d %d %le %le %le %le %le %le %le", &idum1,&idum2,&dum1, &dum2,&dum3,&dum4,&dum5,&dum6,&dum7);
        cell_i[i]=idum1;
		cell_j[i]=idum2;
		T_e[i]=dum1;
		rho[i]=dum2;
		nnh[i]=dum3;
		nne[i]=dum4;
		t[i][0]=dum5;
		t[i][1]=dum6;
		t[i][2]=dum7;	
		v_thermal[i] = pow((2. * KB * T_e[i]/MH),0.5);
//        v_thermal[i]=4.2e5;		
        sigma_e[i] = SIGMA_T * nne[i] / rho[i]; //electron scattering sigma	
        param1[i]=HCLIGHTOVERFOURPI / rho[i] / v_thermal[i];
		
	}
	
	
	if (my_rank==0) printf ("Read physical data - currently using %f Mb\n",g_usedMemory/1e6);
	
	
	return(0);
}
/***************************************************************************/
// Read_cont reads in the continuum flux in the cell. Originally, it was
// set up to read in from a cloudy file - 
/***************************************************************************/


int read_cont()
{  	
  FILE *contf, *output;
  char filename[100];
  char junk[1000];
  double dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9, dum10;
  int idum1,idum2;
  int i,j;
    
	if ((contf = fopen("py_spec_data.dat", "r")) == NULL)
	{
		printf("Cannot open py_spec_data.dat\n");
		abort();
	}

	fscanf(contf, "%*s %d", &idum1);
	model_type=idum1;

    printf ("Model type %i\n",model_type);
	
    
    if (model_type==IONMODE_MATRIX_SPECTRALMODEL)
    {
    	fscanf(contf, "%*s %d", &idum1);
    	nbands=idum1;
    	band_limits=calloc(nbands+1,sizeof(double));	
    	fscanf(contf, "%*s %d", &idum1);
    	if (ncells!=idum1)
    		{
    			printf ("Big problem - number of cells in spec (%i) doesnt equal those in pcon (%i)\n",idum1,ncells);
    			abort();
    		}	
    	for (j=0;j<nbands+1;j++)
    	{
    		fscanf(contf, "%le ",&dum1);
    		band_limits[j]=dum1;
    	}
    	if (my_rank==0) printf ("We have %i J bands running from %e to %e Hz\n",nbands,band_limits[0],band_limits[nbands]);	
    	f1=ARRAY_2D(ncells,nbands,double);
    	f2=ARRAY_2D(ncells,nbands,double);
    	model=ARRAY_2D(ncells,nbands,int);	
    	pl_w=ARRAY_2D(ncells,nbands,double);
    	pl_alpha=ARRAY_2D(ncells,nbands,double);
    	exp_w=ARRAY_2D(ncells,nbands,double);
    	exp_temp=ARRAY_2D(ncells,nbands,double);	
    	mod_fmin=calloc(ncells,sizeof(double));
    	mod_fmax=calloc(ncells,sizeof(double));
    	J=calloc(ncells,sizeof(double));
    	T_max=-1e99; //initialise the maximum temperature
    	T_min=1e99; //initialise the maximum temperature
    	for (i=0;i<ncells;i++)
    	{
    		mod_fmin[i]=1e99;
    		mod_fmax[i]=-1e99;
    		fscanf(contf, " %d", &idum1);
    		fscanf(contf, " %d", &idum2);
    		J[i]=0;
    		for (j=0;j<nbands;j++)
    		{
    			fscanf(contf, "%le %le %d %le %le %le %le",&dum1,&dum2,&idum1,&dum3,&dum4,&dum5,&dum6);
    			f1[i][j]=dum1;
    			f2[i][j]=dum2;
    			model[i][j]=idum1;
    			pl_w[i][j]=dum3;
    			pl_alpha[i][j]=dum4;
    			exp_w[i][j]=dum5;
    			exp_temp[i][j]=dum6;
    			if (model[i][j]!=SPEC_MOD_FAIL && f1[i][j]<mod_fmin[i]) mod_fmin[i]=f1[i][j];
    			if (model[i][j]!=SPEC_MOD_FAIL && f2[i][j]>mod_fmax[i]) mod_fmax[i]=f2[i][j];
    			if (i==0 && my_rank==0) printf ("%e %e %d %e %e %e %e\n",f1[i][j],f2[i][j],model[i][j],pl_w[i][j],pl_alpha[i][j],exp_w[i][j],exp_temp[i][j]);            
    			if (model[i][j]!=SPEC_MOD_FAIL)
    			{
    				J[i]+=int_jnu(i,j);
    			}		
    		}
    //        printf ("%i %e\n",i,J[i]);
    	}
    }
    else if (model_type==IONMODE_MATRIX_BB || model_type == IONMODE_ML93)
    {
    	fscanf(contf, "%*s %d", &idum1);
    	if (ncells!=idum1)
    		{
    			printf ("Big problem - number of cells in spec (%i) doesnt equal those in pcon (%i)\n",idum1,ncells);
    			abort();
    		}        
    	mod_fmin=calloc(ncells,sizeof(double));
    	mod_fmax=calloc(ncells,sizeof(double));
    	W=calloc(ncells,sizeof(double));
    	t_r=calloc(ncells,sizeof(double));
    	J=calloc(ncells,sizeof(double));
        
    	for (i=0;i<ncells;i++)
    	{
    		fscanf(contf, " %d", &idum1);
    		fscanf(contf, " %d", &idum2);
			fscanf(contf, "%le %le ",&dum1,&dum2);
            t_r[i]=dum1;
            W[i]=dum2;
    		mod_fmin[i]=0.0;
    		mod_fmax[i]=100.*KB*t_r[i]/H;
            bb_temp=t_r[i];
            bb_w=W[i];
            J[i]=num_int (bb, mod_fmin[i], mod_fmax[i], 1.e-4);
//            printf ("%i %e %e %e %e\n",i,t_r[i],W[i],mod_fmax[i],J[i]);
            
        }       
    }
	if (my_rank==0) printf ("Read spectral data - currently using %f Mb\n",g_usedMemory/1e6);


	return(0);
}

/***************************************************************************/
//
// Read_line_data reads in all the levels and transitions.
//
/***************************************************************************/

int read_line_data()
{
	FILE *line_list, *atom_models;
	int z,istate,ion_index;
	int i,j,ib;
	int lev_tot;
	int idum1,idum2,idum3;
	double dum1,dum10,dum11,dum2;
	int nlvl,nlevels,nlvl_max;
	int nlines;
    int nline_count;
	
    line_fmin=1e99;
    line_fmax=0.0;
    
    nline_count=0;
	
	if (my_rank==0) printf ("Reading atomic model data\n");	  
	//These lines look thruogh the atomic models file to see how much space we need to assign for models
	if ((atom_models = fopen("atomic_models.txt", "r")) == NULL)  //Try to open atom models data
	{
		printf("Cannot open atomic_models.txt\n");
		abort();
	}
	nlvl_max=-1;	  
	while (fscanf(atom_models, "%d %d %d %lg", &idum1, &idum2, &idum3, &dum1)==4)
	{
		for (i=0;i<nions;i++)
		{
			ion_index=-1;
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_info_nlevels[i]=nlevels=idum3;
				ion_info_xi[i]=dum1;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
			}		
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlevels=idum3;	
		
		for (nlvl = 0; nlvl < nlevels; nlvl++)  //Loop over all the levels for this ion
 		{
			fscanf(atom_models, "%d %lf %lf %d", &idum1, &dum10, &dum11, &idum2); //index - level_energy - multiplicity - unused
 		}
	}
	fclose(atom_models);
	

	lev_energy=ARRAY_2D(nions,nlvl_max,double);
    lev_weight=ARRAY_2D(nions,nlvl_max,double);
    lev_pops=ARRAY_2D(nions,nlvl_max,double); //We will not be filling this array here - but lets set it up while we know nht max number of levels
	
	
	//Now read in the data for real
	
	if ((atom_models = fopen("atomic_models.txt", "r")) == NULL)  //Try to open atom models data
	{
		printf("Cannot open adata.txt.\n");
		abort();
	}
	while (fscanf(atom_models, "%d %d %d %lg", &idum1, &idum2, &idum3, &dum1)==4)
	{
		ion_index=-1;
		for (i=0;i<nions;i++) //Loop over ions to find the matching ion
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_index=i;
			}
			if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
				nlevels=idum3;
			else
			{
				nlevels=ion_info_nlevels[ion_index];
			}
		}
		for (nlvl = 0; nlvl < nlevels; nlvl++)  //Loop over all the levels for this ion
 		{
			fscanf(atom_models, "%d %lf %lf %d", &idum1, &dum10, &dum11, &idum2); //index - level_energy - multiplicity - unused
			if (ion_index !=-1)
			{
				lev_weight[ion_index][nlvl]=dum11;
				lev_energy[ion_index][nlvl]=dum10;	
				
			}
 		}
	}
		

	fclose(atom_models);	
	if (my_rank==0) printf ("Read level data - currently using %f Mb\n",g_usedMemory*1e-6);
		
	
	
	
	//These lines look thruogh the line data file to see how much space we need to assign for models
	if ((line_list = fopen("transitiondata.txt", "r")) == NULL)  //Try to open a file with transition data 
	{
		printf("Cannot open transitiondata.txt\n");
		abort();
	}
	nline_tot=0;
	while (fscanf(line_list, "%d %d %d", &idum1, &idum2, &idum3)==3)
	{
		ion_index=-1;
		for (i=0;i<nions;i++)
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_info_nlines[i]=nlines=idum3;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
				nline_tot=nline_tot+nlines;				
			}		
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlines=idum3;	
		for (nlvl = 0; nlvl < nlines; nlvl++)  //Loop over all the levels for this ion - we dont do anything with this data here
 		{
			fscanf(line_list, "%d %d %d %lf", &idum1, &idum2, &idum3, &dum1); //index - level_energy - multiplicity - unused
 		}
	}
	fclose(line_list);
	
	if (my_rank==0) printf ("There are a total of %i lines linked to ions in our model\n",nline_tot);
	
	//Now we know how many lines we will be reading, we can set up the arrays
		
	trans_ion=calloc(nline_tot,sizeof(int));
	trans_lower=calloc(nline_tot,sizeof(int));
	trans_upper=calloc(nline_tot,sizeof(int));
	trans_iband=calloc(nline_tot,sizeof(int));
		
	trans_freq=calloc(nline_tot,sizeof(double));
	trans_lfreq=calloc(nline_tot,sizeof(double));
	trans_A_ul=calloc(nline_tot,sizeof(double));
	trans_B_ul=calloc(nline_tot,sizeof(double));
	trans_B_lu=calloc(nline_tot,sizeof(double));
	
	g_usedMemory+=4*nline_tot*sizeof(int);
	g_usedMemory+=5*nline_tot*sizeof(double);
	
	
	if ((line_list = fopen("transitiondata.txt", "r")) == NULL)  //Try to open a file with transition data 
	{
		printf("Cannot open transitiondata.txt\n");
		abort();
	}
	nline_tot=0;
	while (fscanf(line_list, "%d %d %d", &idum1, &idum2, &idum3)==3)
	{
		ion_index=-1;
		for (i=0;i<nions;i++)
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				ion_info_nlines[i]=nlines=idum3;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
			}		
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlines=idum3;
		for (nlvl = 0; nlvl < nlines; nlvl++)  //Loop over all the levels for this ion
 		{
			fscanf(line_list, "%d %d %d %lf", &idum1, &idum2, &idum3, &dum1); //index - level_energy - multiplicity - unused
			if (ion_index!=-1) //We only need to add this on if there is a matching ion in our  ion array
			{
				trans_ion[nline_tot]=ion_index;
				if (lev_energy[ion_index][idum3-1] > lev_energy[ion_index][idum2-1])
				{
					trans_lower[nline_tot]=idum2-1; //We subtract 1 because the level array starts at zero
					trans_upper[nline_tot]=idum3-1;
				}
				else //For some reason, the energy levels are swapped round - still use transition but swap the levels
				{
					trans_lower[nline_tot]=idum3-1;
					trans_upper[nline_tot]=idum2-1;
				}
				trans_freq[nline_tot]=(lev_energy[ion_index][trans_upper[nline_tot]]- lev_energy[ion_index][trans_lower[nline_tot]])*EV/H; //Compute the frequency for this transition
				trans_lfreq[nline_tot]=log10(trans_freq[nline_tot]);
	    		trans_A_ul[nline_tot]=dum1; //A_ul
	    		trans_B_ul[nline_tot]=CLIGHTSQUAREDOVERTWOH/pow(trans_freq[nline_tot],3.0) * trans_A_ul[nline_tot]; //B_ul
	    		trans_B_lu[nline_tot]=lev_weight[ion_index][trans_upper[nline_tot]]*trans_B_ul[nline_tot]/
															lev_weight[ion_index][trans_lower[nline_tot]]; //B_luq
				if (trans_freq[nline_tot]<line_fmin) line_fmin=trans_freq[nline_tot];
				if (trans_freq[nline_tot]>line_fmax) line_fmax=trans_freq[nline_tot];
//                trans_iband[nline_tot]=-1; //Set to an error condition                
				for (ib=0;ib<nbands;ib++) //Discover which band this line lies in - saves having to search every time in future.			
				{ 	
//                    trans_iband[nline_tot]=-1; //Set to an error condition
					if (trans_freq[nline_tot]>band_limits[ib] && trans_freq[nline_tot]<band_limits[ib+1])	
                    {
                        trans_iband[nline_tot]=ib;
                        nline_count++;
                    }
				}
				nline_tot++;
			}			
 		}
	}
	fclose(line_list);
	if (my_rank==0) printf ("There are a total of %i lines in modelled bands\n",nline_count);
    
	if (my_rank==0) printf ("Line Fmin=%e Fmax=%e\n",line_fmin,line_fmax);
	if (my_rank==0) printf ("Line Emin=%e Emax=%e\n",line_fmin/(EV/H),line_fmax/(EV/H));
	
	if (my_rank==0) printf ("Read transition data - currently using %f Mb\n",g_usedMemory*1.e-6);	
	
	return(0);
}



/* ********************************************************************* */
char **Array2D (int nx, int ny, size_t dsize)
/*! 
 * Allocate memory for a 2-D array of any basic data type. - stolen from pluto!!!!
 *
 * \param [in] nx    number of elements in the 2nd dimension
 * \param [in] ny    number of elements in the 1st dimension
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char **) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1]
 *          
 *********************************************************************** */
{
  int i;
  char **m;
 
  m    = (char **)calloc ((size_t) nx*sizeof(char *),sizeof(double));
  m[0] = (char *) calloc ((size_t) nx*ny*dsize,sizeof(double));
 
  for (i = 1; i < nx; i++) m[i] = m[(i - 1)] + ny*dsize;
 
  g_usedMemory += nx*ny*dsize;

  return m;
}


double model_jnu(double freq, double lfreq, int icell,int iband,int debug)
{
	int i,j;
	double j_nu;
	
	
    	if (freq>f1[icell][iband] && freq<f2[icell][iband]) //The band could be truncated in this cell
    	{
    		if (model[icell][iband]==SPEC_MOD_FAIL)
    			j_nu=0.0; //No model in this band - should have been caught already but belt and braces
    		else if (model[icell][iband]==SPEC_MOD_PL)
    			j_nu=pow (10, pl_w[icell][iband]+pl_alpha[icell][iband]*lfreq);
    		else if (model[icell][iband]==SPEC_MOD_EXP)
    			j_nu=exp_w[icell][iband] * exp1 ((-1.0 * H * freq) / (KB * exp_temp[icell][iband]));		
    		else 
    		{
    			printf ("Serious problem getting a j_nu for freq %e\n",freq);
    			abort();
    		}				
    	}
    	else
    	{
    		j_nu=0.0; //There is no model - so no flux
    	}

        if (debug==1)
        {
            printf ("Freq %e is in band %i with fmin %e fmax %e model %i\n",freq,iband,f1[icell][iband],f2[icell][iband],model[icell][iband]);
            printf ("exp_w=%e exp_temp=%e exponent=%e\n",exp_w[icell][iband],exp_temp[icell][iband],(-1.0 * H * freq) / (KB * exp_temp[icell][iband]) );
            printf ("exp1=%e exp=%e J_nu=%e\n",exp ((-1.0 * H * freq) / (KB * exp_temp[icell][iband])),exp1 ((-1.0 * H * freq) / (KB * exp_temp[icell][iband])),j_nu);
        }

  
    return(j_nu);
	
		
}


double int_jnu(int icell,int iband)
{
	double a,b,integral,exp_pre;
		
	if (model[icell][iband]==SPEC_MOD_FAIL)
		integral=0.0; //No model in this band - should have been caught already but belt and braces
	else if (model[icell][iband]==SPEC_MOD_PL)
	{
        if (pl_alpha[icell][iband] == -1.0)            //deal with the pathological case
        {
          integral = pow(10.0,pl_w[icell][iband]) * (log(f2[icell][iband])-log(f1[icell][iband]));
        }
        else
        {        
    	    a = pow (10.0, (pl_w[icell][iband] + (log10(f2[icell][iband]) * (pl_alpha[icell][iband] + 1.0))));
    	    b = pow (10.0, (pl_w[icell][iband] + (log10(f1[icell][iband]) * (pl_alpha[icell][iband] + 1.0))));
    	    integral=(a - b) / (pl_alpha[icell][iband] + 1.0);
        }            
	}
	else if (model[icell][iband]==SPEC_MOD_EXP)
	{
	    exp_pre = (-1.0 * H) / (KB * exp_temp[icell][iband]);
		a =1./exp_pre*exp_w[icell][iband] * exp (exp_pre*f2[icell][iband]);
	    b =1./exp_pre*exp_w[icell][iband] * exp (exp_pre*f1[icell][iband]);
		integral=(a - b);		
	}
	else 
	{
		printf ("Serious problem integrating j_nu for band %i cell %i\n",iband,icell);
		abort();
	}
	
	return(integral);
}

double exp1(double x) {
	double ans;
	int indx;
	int test;
	
	if (fabs(x)<0.01)  //At this level - the difference between exp and the approximation is much less than 1pc
	{
		ans=1. + x*(1. + x/2.);
	}
	else if (fabs(x)>500)
	{
		ans=0.0;
	}
	else
	{ 
		ans=exp(x);
	}
		
  
  return ans;
}


double bb(double freq, void *params) {
	double ans;
	int indx;
	int test;
    double bbe;
	
    bbe = exp ((H * freq) / (KB * bb_temp));
    ans = (2. * H * pow (freq, 3.)) / (pow (C, 2));
    ans *= (1 / (bbe - 1));
    ans *= bb_w;  
  return ans;
}


double bb2(double freq, double w, double t_r) {
	double ans;
	int indx;
	int test;
    double bbe;
	
    bbe = exp ((H * freq) / (KB * t_r));
    ans = (2. * H * pow (freq, 3.)) / (pow (C, 2));
    ans *= (1 / (bbe - 1));
    ans *= w;  
  return ans;
}

/**********************************************************/
/**
 * @brief      A wrapper function that carries out numerical integration on a supplied function between limits
 *
 *
 * @param [in] func - a function to integrate, needs to be of the form (double (double,void*))
 * @param [in] a - lower bound
 * @param [in] b - upper bound
 * @param [in] eps - the relative accuracy desired.
 * @return   The integral                          .
 *
 * @details
 * This routine carries out numerical integration of the function func from a to b. It currently
 * uses the gsl implementation of the ROMBERG integration scheme. Initial tests using QAG showed that
 * the type of integrals we do - often with exponential tails - is not well treated. There is an 
 * internal test to ensure that the user isnt asking for an integral where the function is zero everywhere.
 * The gsl routine doesnt handle this well, and takes ages to (sometimes) return zero.
 *
 * ### Notes ###
 *
 **********************************************************/


double
num_int (func, a, b, eps)
     double (*func) (double, void *);
     double a, b;
     double eps;
{
  double result, error, result2;
  double alpha = 0.0;
  void *test = NULL;
  double delta;
  int zflag, i;
  int status = 0;
  int status2 = 0;
  int npoints;
  
  size_t neval;
  gsl_function F;
  F.function = func;
  F.params = &alpha;
  zflag = 1;
  npoints = 1000;
  if (func (a, test) == 0.0 && func (b, test) == 0.0)
  {
    zflag = 0;
    delta = (b - a) / 101;
    for (i = 0; i < 100; i++)
    {
      if (func (a + delta * i, test) != 0)
        zflag = 1;
    }
  }
  
  if (zflag == 1)
  {
    gsl_set_error_handler_off ();       //We need to be able to catch and handle gsl errors 

    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    status = gsl_integration_qags (&F, a, b, 0, eps, 1000, w, &result, &error);
    if (status)
    {
      if (status == GSL_EROUND) //The rounding error has been triggered - try a different integrator
      {
        gsl_integration_workspace_free (w);
        gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc (30);
        status2 = gsl_integration_romberg (&F, a, b, 0, eps, &result2, &neval, w);
        gsl_integration_romberg_free (w);
        if (status2)
        {
          printf ("num_init: some kind of error in romberg and qags integration\n");
        }
      }
    }
    else
    {
      gsl_integration_workspace_free (w);
    }
  }
  else
  {
    result = 0.0;
  }
  return (result);
}



