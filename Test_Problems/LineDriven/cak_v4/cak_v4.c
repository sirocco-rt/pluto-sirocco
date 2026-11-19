#include "cak_v4.h"


/* The multiplier includes ALL eligible lines; no UV-only filter. */

#define NPOINTS 32 //This is the number of t (and hence M) bins
#define tmin -15
#define tmax 5

int my_rank;

/* CAK v4: all-line force multipliers over the supplied spectral range. */
#ifdef MPI_ON
static int mpi_active = 0;
#endif
static int *ion_cell_i, *ion_cell_j;
static void cleanup(void);

static void fail(const char *message)
{
    fprintf(stderr, "cak_v4: %s\n", message);
#ifdef MPI_ON
    if (mpi_active) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static void require_scan(int got, int expected, int line)
{
    if (got != expected) {
        fprintf(stderr, "cak_v4: malformed or truncated input near source line %d (read %d of %d fields)\n", line, got, expected);
        fail("Input parsing failed");
    }
}

static void *checked_calloc(size_t count, size_t size)
{
    if (!count || !size || count > SIZE_MAX / size) fail("Invalid or overflowing allocation size");
    void *ptr = calloc(count, size);
    if (!ptr) fail("Memory allocation failed");
    return ptr;
}
#define calloc checked_calloc

static void require_finite_positive(double x, const char *name)
{
    if (!isfinite(x) || x <= 0.0) fail(name);
}

static void require_end(FILE *fp)
{
    int ch;
    while ((ch = fgetc(fp)) != EOF)
        if (!isspace((unsigned char)ch)) fail("Unexpected trailing input data");
    if (ferror(fp)) fail("Input read failed");
}

static void options(int argc, char **argv)
{
    int i;
    for (i=1; i<argc; ++i) {
        if (!strcmp(argv[i], "--version") || !strcmp(argv[i], "--help")) {
            if (my_rank == 0) {
                puts("cak_v4 4.0 - all-line force multipliers");
                puts("Inputs: sirocco_ion_data.dat, sirocco_pcon_data.dat, sirocco_spec_data.dat");
                puts("Electron opacity: fixed SIROCCO memory2609 sigma_T/m_p (matches the PLUTO line force).");
                puts("Output: M_UV_data.dat (historical name; includes all eligible lines), M_UV_data.meta");
            }
#ifdef MPI_ON
            MPI_Finalize();
#endif
            exit(EXIT_SUCCESS);
        } else fail("Unknown option; use --help");
    }
}

static int same_cell(int row, int i, int j)
{
    return ion_cell_i[row] == i && ion_cell_j[row] == j;
}

static int compare_keys(const void *a, const void *b)
{
    const uint64_t x=*(const uint64_t *)a, y=*(const uint64_t *)b;
    return (x>y)-(x<y);
}

static void check_cell_ids(void)
{
    uint64_t *keys=calloc((size_t)ncells,sizeof(*keys));
    for (int i=0;i<ncells;i++) keys[i]=((uint64_t)(uint32_t)ion_cell_i[i]<<32)|(uint32_t)ion_cell_j[i];
    qsort(keys,(size_t)ncells,sizeof(*keys),compare_keys);
    for (int i=1;i<ncells;i++) if (keys[i]==keys[i-1]) fail("Duplicate cell IDs");
    free(keys);
}


double scale=1.0;

double bb_temp;
double bb_w;

/* Main - top level routine. */
int main(int argc, char** argv)
{

  FILE *output;

  double line_nu;
  double B_lu;
  double kappa_l, delta_doppler, flux_factor, test;
  double stim;
  int upper, lower;

  int icell,itrans,iion,ilev;

//  int my_rank;                  // these two variables are used regardless of parallel mode
  int np_mpi;                   // rank and number of processes, 0 and 1 in non-parallel
  int my_nmin,my_nmax;
  #ifdef MPI_ON
  int np_mpi_global,num_mpi_cells,num_mpi_extra,rank_global;
#endif
  int i,it;
  double *part;
  double t_array[NPOINTS];

#ifdef MPI_ON
  double *M_array_transmit;
  double *M_array_transmit2;
#endif

  double dt;


#ifdef MPI_ON
  MPI_Init (&argc, &argv);
  mpi_active=1;
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif
options(argc, argv);
if (my_rank==0) printf("cak_v4 4.0: all-line calculation; electron opacity=fixed (sigma_T/m_p)\n");
#ifdef MPI_ON
np_mpi_global = np_mpi;
rank_global=my_rank;
#else
(void)np_mpi;
#endif
  g_usedMemory=0.0;



  for (i=0;i<NPOINTS;i++)
  {
	  dt=(tmax-tmin)/(NPOINTS-1.0);
	  t_array[i]=pow(10,tmin+i*dt);
	  if (my_rank==0)
	  {
	  printf ("t[%d] = %e\n",i,t_array[i]);
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

  if (my_rank==0)  printf ("MPI num_mpi_cells=%i num_mpi_extra=%i\n",num_mpi_cells,num_mpi_extra);

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
			if (ion_info_nlevels[iion] == 0) continue; /* Fully stripped/unmodelled ions have no transitions. */
            require_finite_positive(part[iion], "Invalid level partition function");
            part[iion] = 1./part[iion];
			for (ilev=0;ilev<ion_info_nlevels[iion];ilev++) //loop over all the levels for this ion
			{
				lev_pops[iion][ilev]*=part[iion]; //Fraction of the ion in this level (Boltzmann); ion fraction applied in kappa_l
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
					B_lu = trans_B_lu[itrans];
					lower= trans_lower[itrans];
					upper= trans_upper[itrans];
					line_nu= trans_freq[itrans];
					/* Line opacity corrected for stimulated emission:
					     kappa_l ~ B_lu n_l - B_ul n_u
					   With Boltzmann populations and g_l B_lu = g_u B_ul (how B_lu is built
					   in read_line_data) and h nu = E_u - E_l, this is exactly
					     B_lu n_l (1 - exp(-h nu / k T_e)).
					   Evaluating it this way avoids subtracting two nearly equal numbers,
					   so kappa_l can never come out negative through round-off, and the
					   ion fraction (which can be ~1e-300) is multiplied in last so the
					   level populations themselves stay well-scaled. */
					stim = -expm1(-H * line_nu / (KB * T_e[icell]));
					kappa_l = B_lu * lev_pops[iion][lower] * param1[icell];
					kappa_l *= stim;
					kappa_l *= ion_fracs[icell][iion];
					if (!isfinite(kappa_l) || kappa_l < 0.0) /* cannot happen for valid input: report and stop */
					{
						fprintf(stderr,"cak_v4: invalid kappa_l=%e cell=%d (i=%d j=%d) line=%d ion=%d lev %d->%d nu=%e T_e=%e level_frac=%e ion_frac=%e B_lu=%e stim=%e param1=%e\n",
						        kappa_l,icell,cell_i[icell],cell_j[icell],itrans,iion,lower,upper,line_nu,T_e[icell],
						        lev_pops[iion][lower],ion_fracs[icell][iion],B_lu,stim,param1[icell]);
						fail("Invalid line opacity");
					}
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
                            fail("Invalid numerical result; see diagnostic above");
                        }
						if (!isfinite(flux_factor) || flux_factor < 0.0) fail("Invalid spectral weighting");
                        for (it=0;it<NPOINTS;it++)
                        {
	    					if ( (test = kappa_l*t_array[it]/sigma_e[icell]) > 0.0) //Stable expression for every positive optical depth
	    					{
	    						M_UV_array[icell][it] += delta_doppler * flux_factor * (-expm1(-test)) / t_array[it];  //increment the fore multiplier for this t
//	                            by_ion[icell][iion]+=  delta_doppler * flux_factor * (-expm1(-test)) / t_array[it];

	    						if (!isfinite(M_UV_array[icell][it]))
	    						{
	    							printf("Non-finite M1. Abort.\n");
	    							printf("%i %i %i %g %g %g %g\n", icell,cell_i[icell],cell_j[icell],delta_doppler,flux_factor,exp1(-1.*test), t_array[it]);
	    							fail("Invalid numerical result; see diagnostic above");
	    						}
	    					}
	    					else
	    					{
	    						M_UV_array[icell][it] += delta_doppler * flux_factor * kappa_l / sigma_e[icell];
//	                            by_ion[icell][iion]+=  delta_doppler * flux_factor * kappa_l / sigma_e[icell];

	    						if (!isfinite(M_UV_array[icell][it]))
	    						{
	    							printf("Non-finite M2. Abort.\n");
	    							printf("cell=%i ddop=%g ff=%g J=%e kappa_l=%g line_nu=%g\n",icell, delta_doppler,flux_factor,J[icell],kappa_l,line_nu);
	                                if (model_type==IONMODE_MATRIX_SPECTRALMODEL) model_jnu(line_nu, trans_lfreq[itrans],icell,trans_iband[itrans],1);
	    							fail("Invalid numerical result; see diagnostic above");
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

    if ((output = fopen("M_UV_data.dat.tmp", "w")) == NULL)
      {
	printf("Cannot open M__UV_data.dat\n");
	fail("Invalid input; see diagnostic above");
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
		fprintf(output, "%i %i %e %e ", cell_i[icell], cell_j[icell],cell_r[icell],cell_t[icell]);
		for (it=0;it<NPOINTS;it++)
		{
			fprintf (output,"%e ",M_UV_array[icell][it]);
		}
		fprintf(output,"\n");
	}
    if (ferror(output) || fclose(output)!=0) fail("Output write failed");
    if (rename("M_UV_data.dat.tmp","M_UV_data.dat")!=0) fail("Cannot install output table");
    output=fopen("M_UV_data.meta","w");
    if (!output) fail("Cannot write output metadata");
    fprintf(output,"version=cak_v4 4.0\nfrequency_scope=all eligible lines in supplied radiation field\nelectron_opacity=fixed: SIROCCO memory2609 sigma_T/m_p\ninput_density=consistent volume-averaged CMF absolute densities\n");
    if (ferror(output) || fclose(output)!=0) fail("Metadata write failed");

    printf ("output the file\n");




#ifdef MPI_ON
  }
#endif


  free(part);
  cleanup();

#ifdef MPI_ON
  MPI_Finalize ();
  mpi_active=0;
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
  int i,j;
  FILE *abund,*pcon;
  char junk[LINELENGTH];
  double dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8,dum9;
  int idum1,idum2;

/* A loop over the elements - NB at present we need *all* the elements in place - in files which comtain the number of the element as the filename*/


  if ((abund=fopen("sirocco_ion_data.dat","r")) ==NULL)
  {
	  printf ("Cannot open sirocco_ion_data.dat\n");
	  fail("Invalid input; see diagnostic above");
  }

	require_scan(fscanf(abund, "%*s %d", &idum1), 1, __LINE__);
	nions=idum1;
    if (nions <= 0 || nions > 100000) fail("Invalid ion count");
	if (my_rank==0) printf ("We have %i ions\n",nions);

  if (!fgets(junk, sizeof(junk), abund)) fail("Missing ion header terminator");


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
         	fail("Invalid input; see diagnostic above");
        }
		sscanf (junk, "%*s %*d %*s %d %d", &idum1, &idum2);
		if (idum1 <= 0 || idum2 <= 0 || idum2 > idum1+1) fail("Invalid ion identity");
        for (int k=0;k<i;k++) if (ion_info_z[k]==idum1 && ion_info_state[k]==idum2) fail("Duplicate ion identity");
        ion_info_z[i]=idum1;
		ion_info_state[i]=idum2;
	}


	require_scan(fscanf(abund, "%*s %d", &idum1), 1, __LINE__);
	ncells=idum1;
    if (ncells <= 0 || ncells > INT_MAX/NPOINTS) fail("Invalid cell count or MPI reduction count overflow");
    ion_cell_i=calloc((size_t)ncells,sizeof(int));
    ion_cell_j=calloc((size_t)ncells,sizeof(int));
	if (my_rank==0) printf ("We have %i cells\n",ncells);


	ion_fracs=ARRAY_2D(ncells,nions,double); //Set up the ion fracs array


	for (i=0;i<ncells;i++)
	{
		require_scan(fscanf(abund, " %d", &idum1), 1, __LINE__);
		require_scan(fscanf(abund, " %d", &idum2), 1, __LINE__);
        if (idum1 < 0 || idum2 < 0) fail("Negative cell index");
        ion_cell_i[i]=idum1; ion_cell_j[i]=idum2;
		for (j=0;j<nions;j++)
		{
			require_scan(fscanf(abund, "%le",&dum1), 1, __LINE__);
			if (!isfinite(dum1) || dum1 < 0.0) fail("Invalid absolute ion density");
            ion_fracs[i][j]=dum1;
		}
	}
    require_end(abund);
    fclose(abund);
    check_cell_ids();  //Close the file

	if (my_rank==0) printf ("Read ion data - currently using %f Mb\n",g_usedMemory/1e6);



    if ((pcon=fopen("sirocco_pcon_data.dat","r")) ==NULL)
    {
  	  printf ("Cannot open sirocco_pcon_data.dat\n");
  	  fail("Invalid input; see diagnostic above");
    }

	if (fgets (junk, LINELENGTH, pcon) == NULL)
    {
     	printf ("Error reading first line of sirocco_pcon_data.dat\n");
     	fail("Invalid input; see diagnostic above");
    }

	require_scan(sscanf(junk, "%*s %d", &idum1), 1, __LINE__);


	if (idum1!=ncells)
	{
		printf ("py_pcon_data has different number of cells - aborting");
		fail("Invalid input; see diagnostic above");
	}

	cell_j=calloc(ncells,sizeof(int));
	cell_i=calloc(ncells,sizeof(int));
	cell_r=calloc(ncells,sizeof(double));
	cell_t=calloc(ncells,sizeof(double));
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
         	fail("Invalid input; see diagnostic above");
        }
		sscanf (junk, "%d %d %le %le %le %le %le %le %le %le %le", &idum1,&idum2,&dum1, &dum2,&dum3,&dum4,&dum5,&dum6,&dum7,&dum8, &dum9);
        if (!same_cell(i,idum1,idum2)) fail("Ion and physical cell IDs/order differ");
        cell_i[i]=idum1;
		cell_j[i]=idum2;
		cell_r[i]=dum1;
		cell_t[i]=dum2;
		T_e[i]=dum3;
		rho[i]=dum4;
		nnh[i]=dum5;
		nne[i]=dum6;
		t[i][0]=dum7;
		t[i][1]=dum8;
		t[i][2]=dum9;
		require_finite_positive(T_e[i], "Electron temperature must be positive and finite");
        require_finite_positive(rho[i], "Mass density must be positive and finite (check partial/empty cells)");
        require_finite_positive(nne[i], "Electron density must be positive and finite");
        if (!isfinite(nnh[i]) || nnh[i]<0 || !isfinite(cell_r[i]) || !isfinite(cell_t[i])) fail("Invalid physical cell data");
        for (int band=0;band<3;band++) if (!isfinite(t[i][band]) || t[i][band]<0) fail("Invalid supplied optical-depth parameter");
        v_thermal[i] = pow((2. * KB * T_e[i]/MH),0.5);
//        v_thermal[i]=4.2e5;
        sigma_e[i] = SIROCCO_FIXED_SIGMA_E; //fixed electron scattering opacity sigma_T/m_p (matches PLUTO)
        param1[i]=HCLIGHTOVERFOURPI / rho[i] / v_thermal[i];
        require_finite_positive(sigma_e[i], "Invalid electron opacity");
        require_finite_positive(param1[i], "Invalid opacity normalization");

	}


	require_end(pcon);
    fclose(pcon);
	if (my_rank==0) printf ("Read physical data - currently using %f Mb\n",g_usedMemory/1e6);


	return(0);
}
/***************************************************************************/
// Read_cont reads in the continuum flux in the cell. Originally, it was
// set up to read in from a cloudy file -
/***************************************************************************/


int read_cont()
{
  FILE *contf;
  double dum1, dum2, dum3, dum4, dum5, dum6;
  int idum1,idum2;
  int i,j;
  long n_exp_bad = 0, n_exp_negT = 0;

	if ((contf = fopen("sirocco_spec_data.dat", "r")) == NULL)
	{
		printf("Cannot open sirocco_spec_data.dat\n");
		fail("Invalid input; see diagnostic above");
	}

	require_scan(fscanf(contf, "%*s %d", &idum1), 1, __LINE__);
	model_type=idum1;
    if (model_type!=IONMODE_MATRIX_SPECTRALMODEL && model_type!=IONMODE_MATRIX_BB && model_type!=IONMODE_ML93) fail("Unsupported radiation model");

    printf ("Model type %i\n",model_type);


    if (model_type==IONMODE_MATRIX_SPECTRALMODEL)
    {
    	require_scan(fscanf(contf, "%*s %d", &idum1), 1, __LINE__);
    	nbands=idum1;
        if (nbands<=0 || nbands>100000) fail("Invalid spectral band count");
    	band_limits=calloc(nbands+1,sizeof(double));
    	require_scan(fscanf(contf, "%*s %d", &idum1), 1, __LINE__);
    	if (ncells!=idum1)
    		{
    			printf ("Big problem - number of cells in spec (%i) doesnt equal those in pcon (%i)\n",idum1,ncells);
    			fail("Invalid input; see diagnostic above");
    		}
    	for (j=0;j<nbands+1;j++)
    	{
    		require_scan(fscanf(contf, "%le ",&dum1), 1, __LINE__);
    		require_finite_positive(dum1,"Invalid spectral band boundary");
            if (j && dum1<=band_limits[j-1]) fail("Spectral boundaries must increase");
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
    		require_scan(fscanf(contf, " %d", &idum1), 1, __LINE__);
    		require_scan(fscanf(contf, " %d", &idum2), 1, __LINE__);
            if (!same_cell(i,idum1,idum2)) fail("Ion and spectral cell IDs/order differ");
    		J[i]=0;
    		for (j=0;j<nbands;j++)
    		{
    			require_scan(fscanf(contf, "%le %le %d %le %le %le %le",&dum1,&dum2,&idum1,&dum3,&dum4,&dum5,&dum6), 7, __LINE__);
    			f1[i][j]=dum1;
    			f2[i][j]=dum2;
    			model[i][j]=idum1;
    			pl_w[i][j]=dum3;
    			pl_alpha[i][j]=dum4;
    			exp_w[i][j]=dum5;
    			exp_temp[i][j]=dum6;
    			if (model[i][j]!=SPEC_MOD_FAIL) {
                    if (model[i][j]!=SPEC_MOD_PL && model[i][j]!=SPEC_MOD_EXP) fail("Unsupported spectral fit type");
                    require_finite_positive(dum1,"Invalid spectral minimum");
                    require_finite_positive(dum2,"Invalid spectral maximum");
                    if (dum2<=dum1 || dum1<band_limits[j]*(1-2e-6) || dum2>band_limits[j+1]*(1+2e-6)) fail("Cell spectral fit lies outside its global band");
                    if (model[i][j]==SPEC_MOD_PL && (!isfinite(dum3)||!isfinite(dum4))) fail("Invalid power-law fit");
                    if (model[i][j]==SPEC_MOD_EXP) {
                        /* j_nu = w*exp(-h nu / k T). A negative finite T is a valid
                           (rising) exponential fit, so keep it. Only a zero/NaN/inf T
                           or a zero/NaN/negative weight is unusable: treat that band as
                           "no model" (zero contribution) and report it, instead of
                           aborting the whole MPI job. */
                        int bad_w = (!isfinite(dum5) || dum5 < 0.0);
                        int bad_t = (!isfinite(dum6) || dum6 == 0.0);
                        if (dum6 < 0.0 && isfinite(dum6)) n_exp_negT++;
                        if (bad_w || bad_t || dum5 == 0.0) {
                            if (my_rank==0 && n_exp_bad < 20)
                                fprintf(stderr, "cak_v4: WARNING cell %d band %d: exponential fit w=%e T=%e unusable -> band set to no-model\n",
                                        i, j, dum5, dum6);
                            n_exp_bad++;
                            model[i][j]=SPEC_MOD_FAIL;
                        }
                    }
                }
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
    	require_scan(fscanf(contf, "%*s %d", &idum1), 1, __LINE__);
    	if (ncells!=idum1)
    		{
    			printf ("Big problem - number of cells in spec (%i) doesnt equal those in pcon (%i)\n",idum1,ncells);
    			fail("Invalid input; see diagnostic above");
    		}
    	mod_fmin=calloc(ncells,sizeof(double));
    	mod_fmax=calloc(ncells,sizeof(double));
    	W=calloc(ncells,sizeof(double));
    	t_r=calloc(ncells,sizeof(double));
    	J=calloc(ncells,sizeof(double));

    	for (i=0;i<ncells;i++)
    	{
    		require_scan(fscanf(contf, " %d", &idum1), 1, __LINE__);
    		require_scan(fscanf(contf, " %d", &idum2), 1, __LINE__);
            if (!same_cell(i,idum1,idum2)) fail("Ion and spectral cell IDs/order differ");
			require_scan(fscanf(contf, "%le %le ",&dum1,&dum2), 2, __LINE__);
            require_finite_positive(dum1,"Invalid radiation temperature");
            if (!isfinite(dum2) || dum2<0) fail("Invalid blackbody dilution");
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
	require_end(contf);
    fclose(contf);
    for (i=0;i<ncells;i++) if (!isfinite(J[i]) || J[i]<0.0) fail("Invalid integrated radiation intensity");
    if (my_rank==0 && (n_exp_bad || n_exp_negT))
        printf ("Spectral model check: %ld exponential band(s) unusable (set to no-model), %ld with negative T (kept)\n", n_exp_bad, n_exp_negT);
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
	int ion_index = -1;
	int i,ib;
		int idum1,idum2,idum3;
	double dum1,dum10,dum11;
	int nlvl,nlevels=0,nlvl_max;
	int nlines;
    int nline_count, header_count;
    int *seen_models=calloc((size_t)nions,sizeof(int));
    int *seen_lines=calloc((size_t)nions,sizeof(int));

    line_fmin=1e99;
    line_fmax=0.0;

    nline_count=0;

	if (my_rank==0) printf ("Reading atomic model data\n");
	//These lines look thruogh the atomic models file to see how much space we need to assign for models
	if ((atom_models = fopen("atomic_models.txt", "r")) == NULL)  //Try to open atom models data
	{
		printf("Cannot open atomic_models.txt\n");
		fail("Invalid input; see diagnostic above");
	}
	nlvl_max=-1;
	while ((header_count=fscanf(atom_models, "%d %d %d %lg", &idum1, &idum2, &idum3, &dum1))==4)
	{
        if (idum1<=0 || idum2<=0 || idum3<0 || idum3>1000000) fail("Invalid atomic model header");
		for (i=0;i<nions;i++)
		{
			ion_index=-1;
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				if (idum3<0 || idum3>1000000) fail("Invalid atomic level count");
                if (seen_models[i]++) fail("Duplicate atomic model ion");
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
			require_scan(fscanf(atom_models, "%d %lf %lf %d", &idum1, &dum10, &dum11, &idum2), 4, __LINE__); //index - level_energy - multiplicity - unused
 		}
	}
	if (header_count!=EOF) fail("Malformed atomic model header");
    fclose(atom_models);


	if (nlvl_max<=0) fail("No atomic models match the supplied ions");
    lev_energy=ARRAY_2D(nions,nlvl_max,double);
    lev_weight=ARRAY_2D(nions,nlvl_max,double);
    lev_pops=ARRAY_2D(nions,nlvl_max,double); //We will not be filling this array here - but lets set it up while we know nht max number of levels


	//Now read in the data for real

	if ((atom_models = fopen("atomic_models.txt", "r")) == NULL)  //Try to open atom models data
	{
		printf("Cannot open adata.txt.\n");
		fail("Invalid input; see diagnostic above");
	}
	while ((header_count=fscanf(atom_models, "%d %d %d %lg", &idum1, &idum2, &idum3, &dum1))==4)
	{
        if (idum1<=0 || idum2<=0 || idum3<0 || idum3>1000000) fail("Invalid atomic model header");
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
			require_scan(fscanf(atom_models, "%d %lf %lf %d", &idum1, &dum10, &dum11, &idum2), 4, __LINE__); //index - level_energy - multiplicity - unused
			if (ion_index !=-1)
			{
				if (idum1!=nlvl+1 || !isfinite(dum10) || dum10<0 || !isfinite(dum11) || dum11<=0) fail("Invalid atomic level record");
                lev_weight[ion_index][nlvl]=dum11;
				lev_energy[ion_index][nlvl]=dum10;

			}
 		}
	}


	if (header_count!=EOF) fail("Malformed atomic model header");
    fclose(atom_models);
	if (my_rank==0) printf ("Read level data - currently using %f Mb\n",g_usedMemory*1e-6);




	//These lines look thruogh the line data file to see how much space we need to assign for models
	if ((line_list = fopen("transitiondata.txt", "r")) == NULL)  //Try to open a file with transition data
	{
		printf("Cannot open transitiondata.txt\n");
		fail("Invalid input; see diagnostic above");
	}
	nline_tot=0;
	while ((header_count=fscanf(line_list, "%d %d %d", &idum1, &idum2, &idum3))==3)
	{
        if (idum1<=0 || idum2<=0 || idum3<0 || idum3>INT_MAX-nline_tot) fail("Invalid transition header");
		ion_index=-1;
		for (i=0;i<nions;i++)
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				if (idum3<0 || idum3>INT_MAX-nline_tot) fail("Invalid transition count");
                if (seen_lines[i]++) fail("Duplicate transition ion block");
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
			require_scan(fscanf(line_list, "%d %d %d %lf", &idum1, &idum2, &idum3, &dum1), 4, __LINE__); //index - level_energy - multiplicity - unused
 		}
	}
	if (header_count!=EOF) fail("Malformed transition header");
    fclose(line_list);

	if (my_rank==0) printf ("There are a total of %i lines linked to ions in our model\n",nline_tot);

	//Now we know how many lines we will be reading, we can set up the arrays

	if (nline_tot<=0) fail("No transitions match the supplied ions");
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
		fail("Invalid input; see diagnostic above");
	}
	nline_tot=0;
	while ((header_count=fscanf(line_list, "%d %d %d", &idum1, &idum2, &idum3))==3)
	{
        if (idum1<=0 || idum2<=0 || idum3<0 || idum3>INT_MAX-nline_tot) fail("Invalid transition header");
		ion_index=-1;
		for (i=0;i<nions;i++)
		{
			if (idum1==ion_info_z[i] && idum2==ion_info_state[i]) //Find the correct ion record
			{
				if (idum3<0 || idum3>INT_MAX-nline_tot) fail("Invalid transition count");
                ion_info_nlines[i]=nlines=idum3;
				ion_index=i;
				if (nlevels>nlvl_max) nlvl_max=nlevels;
			}
		}
		if (ion_index==-1) //We didnt match an ion in our data - we still need to skim over tho
			nlines=idum3;
		for (nlvl = 0; nlvl < nlines; nlvl++)  //Loop over all the levels for this ion
 		{
			require_scan(fscanf(line_list, "%d %d %d %lf", &idum1, &idum2, &idum3, &dum1), 4, __LINE__); //index - level_energy - multiplicity - unused
			if (ion_index!=-1) //We only need to add this on if there is a matching ion in our  ion array
			{
				if (idum2<1 || idum3<1 || idum2>ion_info_nlevels[ion_index] || idum3>ion_info_nlevels[ion_index] || !isfinite(dum1) || dum1<0) fail("Invalid transition level index or Einstein A");
                if (lev_energy[ion_index][idum2-1]==lev_energy[ion_index][idum3-1]) fail("Zero-energy transition");
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
trans_iband[nline_tot]=(model_type==IONMODE_MATRIX_SPECTRALMODEL)?-1:0;
				for (ib=0;ib<nbands;ib++) //Discover which band this line lies in - saves having to search every time in future.
				{
//                    trans_iband[nline_tot]=-1; //Set to an error condition
					if (trans_freq[nline_tot]>=band_limits[ib] && (trans_freq[nline_tot]<band_limits[ib+1] || (ib==nbands-1 && trans_freq[nline_tot]==band_limits[ib+1])))
                    {
                        trans_iband[nline_tot]=ib;
                        nline_count++;
                    }
				}
				nline_tot++;
			}
 		}
	}
	if (header_count!=EOF) fail("Malformed transition header");
    fclose(line_list);
	if (my_rank==0) printf ("There are a total of %i lines in modelled bands\n",nline_count);

	if (my_rank==0) printf ("Line Fmin=%e Fmax=%e\n",line_fmin,line_fmax);
	if (my_rank==0) printf ("Line Emin=%e Emax=%e\n",line_fmin/(EV/H),line_fmax/(EV/H));

	free(seen_models); free(seen_lines);
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

  if (nx<=0 || ny<=0 || (size_t)nx>SIZE_MAX/(size_t)ny) fail("Invalid 2-D allocation");
  m    = (char **)calloc ((size_t) nx,sizeof(char *));
  m[0] = (char *) calloc ((size_t) nx*(size_t)ny,dsize);

  for (i = 1; i < nx; i++) m[i] = m[(i - 1)] + ny*dsize;

  g_usedMemory += nx*ny*dsize;

  return m;
}


double model_jnu(double freq, double lfreq, int icell,int iband,int debug)
{
	double j_nu;


    	if (freq>=f1[icell][iband] && freq<=f2[icell][iband]) //The band could be truncated in this cell
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
    			fail("Invalid input; see diagnostic above");
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
    long double a=f1[icell][iband], b=f2[icell][iband], result=0;
    if (model[icell][iband]==SPEC_MOD_FAIL) return 0.0;
    if (model[icell][iband]==SPEC_MOD_PL) {
        long double p=(long double)pl_alpha[icell][iband]+1.0L;
        long double log_ratio=logl(b/a);
        long double log_norm=logl(10.0L)*(long double)pl_w[icell][iband];
        if (p==0) result=expl(log_norm)*log_ratio;
        else if (p>0) result=expl(log_norm+p*logl(b))*(-expm1l(-p*log_ratio))/p;
        else result=expl(log_norm+p*logl(a))*expm1l(p*log_ratio)/p;
    } else if (model[icell][iband]==SPEC_MOD_EXP) {
        long double q=(long double)H/((long double)KB*exp_temp[icell][iband]);
        result=(long double)exp_w[icell][iband]*expl(-q*a)*(-expm1l(-q*(b-a)))/q;
    } else fail("Unsupported spectral fit in integral");
    if (!isfinite((double)result) || result<0) fail("Invalid spectral integral");
    return (double)result;
}

double exp1(double x) { return exp(x); }


double bb(double freq, void *params) {
    (void)params;
	double ans;
    double bbe;

    if (freq<=0 || bb_w==0) return 0.0;
    bbe = expm1 ((H * freq) / (KB * bb_temp));
    ans = (2. * H * pow (freq, 3.)) / (pow (C, 2));
    ans *= (1 / bbe);
    ans *= bb_w;
  return ans;
}


double bb2(double freq, double w, double t_r) {
	double ans;
    double bbe;

    if (freq<=0 || w==0) return 0.0;
    bbe = expm1 ((H * freq) / (KB * t_r));
    ans = (2. * H * pow (freq, 3.)) / (pow (C, 2));
    ans *= (1 / bbe);
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
num_int(double (*func)(double,void *), double a, double b, double eps)
{
    gsl_function F;
    double result=0.0,error=0.0,alpha=0.0;
    size_t neval;
    int status;
    F.function=func; F.params=&alpha;
    if (!isfinite(a)||!isfinite(b)||b<=a||eps<=0) fail("Invalid integration interval");
    gsl_set_error_handler_off();
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    if (!w) fail("Cannot allocate integration workspace");
    status=gsl_integration_qags(&F,a,b,0,eps,1000,w,&result,&error);
    gsl_integration_workspace_free(w);
    if (status==GSL_EROUND) {
        gsl_integration_romberg_workspace *rw=gsl_integration_romberg_alloc(20);
        if (!rw) fail("Cannot allocate fallback integration workspace");
        status=gsl_integration_romberg(&F,a,b,0,eps,&result,&neval,rw);
        gsl_integration_romberg_free(rw);
    }
    if (status!=GSL_SUCCESS || !isfinite(result)) fail("Radiation integration failed");
    return result;
}

static void free2d(void *array)
{
    void **p=array;
    if (p) {free(p[0]); free(p);}
}

static void cleanup(void)
{
    free(ion_info_z); free(ion_info_state); free(ion_info_nlevels); free(ion_info_nlines); free(ion_info_xi);
    free(ion_cell_i); free(ion_cell_j); free(cell_i); free(cell_j); free(cell_r); free(cell_t);
    free(T_e); free(nnh); free(nne); free(rho); free(v_thermal); free(sigma_e); free(param1);
    free2d(ion_fracs); free2d(by_ion); free2d(t); free2d(M_UV_array);
    free2d(lev_energy); free2d(lev_weight); free2d(lev_pops);
    free2d(f1); free2d(f2); free2d(model); free2d(pl_w); free2d(pl_alpha); free2d(exp_w); free2d(exp_temp);
    free(mod_fmin); free(mod_fmax); free(band_limits); free(J); free(W); free(t_r);
    free(trans_ion); free(trans_lower); free(trans_upper); free(trans_iband);
    free(trans_freq); free(trans_lfreq); free(trans_A_ul); free(trans_B_ul); free(trans_B_lu);
}
