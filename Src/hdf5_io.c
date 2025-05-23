/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief HDF5 I/O main driver for static grid implementation.

  This file provides I/O functionality for HDF5 data format.
  The output suffix is ".h5".

  The WriteHDF5() function allows to write data in serial or parallel
  mode in either single or double precision.

  The ReadHDF5() function allows to read double precision data 
  in serial or parrallel mode.
  
  \note 
  By turning "MPI_POSIX" to "YES", HDF5 uses another parallel IO driver
  called MPI POSIX which is a "combination" MPI-2 and posix I/O driver.
  It uses MPI for coordinating the actions of several processes and
  posix I/O calls to do the actual I/O to the disk.
  There is no collective I/O mode with this driver. This will almost
  certainly not work correctly for files accessed on distributed parallel
  systems with the file located on a non-parallel filesystem.
  On some systems, Using MPI POSIX driver may perform better than
  using MPI-IO driver with independent IO mode.\n
  For more info take a look at
  http://www.hdfgroup.org/HDF5/PHDF5/parallelhdf5hints.pdf

  \authors C. Zanni (zanni@oato.inaf.it)\n
           A. Mignone (andrea.mignone@unito.it)
           G. Musicanisi (g.muscianisi@cineca.it)\n
  \date   Jul 21, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#define H5_USE_16_API
#include "hdf5.h"

#ifndef PARALLEL
 #define MPI_POSIX YES
#else
 #define MPI_POSIX NO
#endif

/* ********************************************************************* */
void WriteHDF5 (Output *output, Grid *grid)
/*!
 * Write data to disk using hdf5 format in single or double
 * precision.
 *
 * \param [in] output the output structure associated with HDF5 format
 * \param [in] grid   a pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  hid_t dataspace, memspace, dataset;
  hid_t strspace, stratt, string_type;
  hid_t tspace, tattr;
  hid_t file_identifier, group, timestep;
  hid_t file_access = 0;
 #if MPI_POSIX == NO
  hid_t plist_id_mpiio = 0; /* for collective MPI I/O */
 #endif
  hid_t err;

  hsize_t dimstr;
  hsize_t dimens[DIMENSIONS];
  hsize_t start[DIMENSIONS];
  hsize_t stride[DIMENSIONS];
  hsize_t count[DIMENSIONS];

  time_t tbeg, tend;
  static float ****node_coords, ****cell_coords;
  char filename[512], filenamexmf[512], tstepname[32];
  char *coords = "/cell_coords/X /cell_coords/Y /cell_coords/Z ";
  char *cname[] = {"X", "Y", "Z"};
  char xmfext[8];
  int  rank, nd, nr, nv, ns, nc, ngh, ii, jj, kk;
  int n1p, n2p, n3p, nprec;
  FILE *fxmf;

/* ----------------------------------------------------------------
                 compute coordinates just once
   ---------------------------------------------------------------- */

  if (node_coords == NULL) {
    double x1, x2, x3;

    n1p = NX1 + (grid->rbound[IDIR] != 0);
    n2p = NX2 + (grid->rbound[JDIR] != 0);
    n3p = NX3 + (grid->rbound[KDIR] != 0);

    node_coords = ARRAY_4D(3, n3p, n2p, n1p, float);
    cell_coords = ARRAY_4D(3, NX3, NX2, NX1, float);
    
    for (kk = 0; kk < n3p; kk++) {  x3 = grid->xl[KDIR][KBEG+kk];
    for (jj = 0; jj < n2p; jj++) {  x2 = grid->xl[JDIR][JBEG+jj];
    for (ii = 0; ii < n1p; ii++) {  x1 = grid->xl[IDIR][IBEG+ii];

      node_coords[JDIR][kk][jj][ii] = 0.0;
      node_coords[KDIR][kk][jj][ii] = 0.0;
     
      #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
      DIM_EXPAND(node_coords[IDIR][kk][jj][ii] = (float)x1;  ,
               node_coords[JDIR][kk][jj][ii] = (float)x2;  ,
               node_coords[KDIR][kk][jj][ii] = (float)x3;)
      #elif GEOMETRY == POLAR
      DIM_EXPAND(node_coords[IDIR][kk][jj][ii] = (float)(x1*cos(x2)); ,
               node_coords[JDIR][kk][jj][ii] = (float)(x1*sin(x2)); ,
               node_coords[KDIR][kk][jj][ii] = (float)(x3);)
      #elif GEOMETRY == SPHERICAL
        #if DIMENSIONS <= 2
        DIM_EXPAND(node_coords[IDIR][kk][jj][ii] = (float)(x1*sin(x2)); ,
                 node_coords[JDIR][kk][jj][ii] = (float)(x1*cos(x2)); ,
                 node_coords[KDIR][kk][jj][ii] = 0.0;)
        #else
        DIM_EXPAND(node_coords[IDIR][kk][jj][ii] = (float)(x1*sin(x2)*cos(x3)); ,
                 node_coords[JDIR][kk][jj][ii] = (float)(x1*sin(x2)*sin(x3)); ,
                 node_coords[KDIR][kk][jj][ii] = (float)(x1*cos(x2));)
        #endif
      #else
      printLog ("! HDF5_IO: Unknown geometry\n");
      QUIT_PLUTO(1);
      #endif
    }}}

    for (kk = 0; kk < NX3; kk++) {  x3 = grid->x[KDIR][KBEG+kk];
    for (jj = 0; jj < NX2; jj++) {  x2 = grid->x[JDIR][JBEG+jj];
    for (ii = 0; ii < NX1; ii++) {  x1 = grid->x[IDIR][IBEG+ii];

      cell_coords[JDIR][kk][jj][ii] = 0.0;
      cell_coords[KDIR][kk][jj][ii] = 0.0;

      #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
       DIM_EXPAND(cell_coords[IDIR][kk][jj][ii] = (float)x1;  ,
                cell_coords[JDIR][kk][jj][ii] = (float)x2;  ,
                cell_coords[KDIR][kk][jj][ii] = (float)x3;)
      #elif GEOMETRY == POLAR
       DIM_EXPAND(cell_coords[IDIR][kk][jj][ii] = (float)(x1*cos(x2)); ,
                cell_coords[JDIR][kk][jj][ii] = (float)(x1*sin(x2)); ,
                cell_coords[KDIR][kk][jj][ii] = (float)(x3);)
      #elif GEOMETRY == SPHERICAL
       #if DIMENSIONS <= 2
        DIM_EXPAND(cell_coords[IDIR][kk][jj][ii] = (float)(x1*sin(x2)); ,
                 cell_coords[JDIR][kk][jj][ii] = (float)(x1*cos(x2)); ,
                 cell_coords[KDIR][kk][jj][ii] = 0.0;)
       #else
        DIM_EXPAND(cell_coords[IDIR][kk][jj][ii] = (float)(x1*sin(x2)*cos(x3)); ,
                 cell_coords[JDIR][kk][jj][ii] = (float)(x1*sin(x2)*sin(x3)); ,
                 cell_coords[KDIR][kk][jj][ii] = (float)(x1*cos(x2));)
       #endif
      #else
       printLog ("! HDF5_IO: Unknown geometry\n");
       QUIT_PLUTO(1);
      #endif
    }}}
  } 

/* --------------------------------------------------------------
   ! Important: data is written in reverse order (Z-Y-X),
     and we will often use nr = DIMENSIONS-nd-1
   -------------------------------------------------------------- */

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (prank == 0)time(&tbeg);
#endif

  sprintf (filename, "%s/data.%04d.%s", output->dir,output->nfile, output->ext);

  rank   = DIMENSIONS;
  dimstr = 3;

#ifdef PARALLEL
  file_access = H5Pcreate(H5P_FILE_ACCESS);
  #if MPI_POSIX == YES
  H5Pset_fapl_mpiposix(file_access, MPI_COMM_WORLD, 1); 
  #else
  H5Pset_fapl_mpio(file_access,  MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  file_identifier = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, file_access);
  H5Pclose(file_access);
#else
  file_identifier = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif

  sprintf (tstepname, "Timestep_%d", output->nfile);
  timestep = H5Gcreate(file_identifier, tstepname, 0);

  tspace  = H5Screate(H5S_SCALAR);
  tattr   = H5Acreate(timestep, "Time", H5T_NATIVE_DOUBLE, tspace, H5P_DEFAULT);
  err = H5Awrite(tattr, H5T_NATIVE_DOUBLE, &g_time);
  H5Aclose(tattr);
  H5Sclose(tspace);

  group = H5Gcreate(timestep, "vars", 0); /* Create group "vars" (cell-centered vars) */

/* Define "coords" attribute of group "vars" */

  strspace = H5Screate_simple(1, &dimstr, NULL);
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen("/cell_coords/X "));
  H5Tset_strpad(string_type,H5T_STR_SPACEPAD);
  stratt = H5Acreate(group,"coords",string_type,strspace, H5P_DEFAULT);
  err = H5Awrite(stratt, string_type, coords);
  H5Aclose(stratt);
  H5Sclose(strspace);

  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    dimens[nd] = grid->np_int_glob[nr];
  }
 
  dataspace = H5Screate_simple(rank, dimens, NULL);

#ifdef PARALLEL 
  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    start[nd]  = grid->beg[nr] - grid->nghost[nr];
    stride[nd] = 1;
    count[nd]  = grid->np_int[nr];
  }
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, NULL);
#endif

  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    dimens[nd] = grid->np_tot[nr];
  }

  memspace = H5Screate_simple(rank,dimens,NULL);

  for (nd = 0; nd < DIMENSIONS; nd++){
    nr = DIMENSIONS-nd-1;
    start[nd]  = grid->nghost[nr];
    stride[nd] = 1;
    count[nd]  = grid->np_int[nr];
  }
  err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET, start, stride, count, NULL);

/* ------------------------------------
      write cell-centered data
   ------------------------------------ */

#if MPI_POSIX == NO
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
#endif

  for (nv = 0; nv < output->nvar; nv++) {
  
  /* -- skip variable if excluded from output or if it is staggered -- */

    if (!output->dump_var[nv] || output->stag_var[nv] != -1) continue;

    if (output->type == DBL_H5_OUTPUT){
      dataset = H5Dcreate(group, output->var_name[nv], H5T_NATIVE_DOUBLE,
                          dataspace, H5P_DEFAULT);
      #if MPI_POSIX == NO
      err = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     plist_id_mpiio, output->V[nv][0][0]);
      #else
      err = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, output->V[nv][0][0]);
      #endif
    }else if (output->type == FLT_H5_OUTPUT){
      void *Vpt;
      Vpt = (void *)(Convert_dbl2flt(output->V[nv],1.0, 0))[0][0];

      dataset = H5Dcreate(group, output->var_name[nv], H5T_NATIVE_FLOAT,
                          dataspace, H5P_DEFAULT);

      #if MPI_POSIX == NO
      err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                     dataspace, plist_id_mpiio, Vpt);
      #else
      err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                     dataspace, H5P_DEFAULT, Vpt);
      #endif
    }
    H5Dclose(dataset);
  }

#if MPI_POSIX == NO
  H5Pclose(plist_id_mpiio);
#endif
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Gclose(group); /* Close group "vars" */

/*  --------------------------------------------------------------------- */
/*! \note Writing of staggered magnetic fields:
    Here we follow the same convention used by ArrayLib, where internal
    points that lie on the grid boundary between processors are owned
    by the \b left processor.
    For this reason, only the leftmost processor writes one additional
    point in the direction of staggering while the other processors
    write the same amount of zones as for cell-centered data.
    Also, keep in mind that the data array Vs starts from -1
    (instead of 0) in the staggered direction.                           */
/*  --------------------------------------------------------------------- */

#ifdef STAGGERED_MHD
  if (output->type == DBL_H5_OUTPUT){
    group = H5Gcreate(timestep, "stag_vars", 0); /* Create group "stag_vars" (staggered vars) */

    #if MPI_POSIX == NO
    plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
    #endif

    for (ns = 0; ns < DIMENSIONS; ns++) {

  /* -- skip variable if excluded from output or if it is cell-centered -- */
 
      if (!output->dump_var[NVAR+ns] || output->stag_var[NVAR+ns] == -1) continue;
      for (nd = 0; nd < DIMENSIONS; nd++) {
        nr = DIMENSIONS-nd-1;
        dimens[nd] = grid->np_int_glob[nr] + (ns == nr);
      }
      dataspace = H5Screate_simple(rank, dimens, NULL);

      #ifdef PARALLEL
      for (nd = 0; nd < DIMENSIONS; nd++) {
        nr = DIMENSIONS-nd-1;
        start[nd]  = grid->beg[nr] - grid->nghost[nr];
        stride[nd] = 1;
        count[nd]  = grid->np_int[nr];
        if (ns == nr){
           if (grid->lbound[ns] != 0) count[nd] += 1;
           else                       start[nd] += 1;
        }
      }
      err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                                start, stride, count, NULL);
      #endif

      for (nd = 0; nd < DIMENSIONS; nd++){
        nr = DIMENSIONS-nd-1;
        dimens[nd] = grid->np_tot[nr] + (ns==nr);
      }
      memspace = H5Screate_simple(rank,dimens,NULL);

      for (nd = 0; nd < DIMENSIONS; nd++){
        nr = DIMENSIONS-nd-1;
        start[nd]  = grid->nghost[nr];
        stride[nd] = 1;
        count[nd]  = grid->np_int[nr];
        if (ns == nr && grid->lbound[ns] != 0) {
          start[nd] -= 1;
          count[nd] += 1;
        }
      }

      err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET, start, stride, count, NULL);

      if (output->type == DBL_H5_OUTPUT){
        dataset = H5Dcreate(group, output->var_name[NVAR+ns], H5T_NATIVE_DOUBLE,
                           dataspace, H5P_DEFAULT);
        #if MPI_POSIX == NO
        err = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                       plist_id_mpiio, output->V[NVAR+ns][0][0]);
        #else
        err = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                       H5P_DEFAULT, output->V[NVAR+ns][0][0]);
        #endif
      }else if (output->type == FLT_H5_OUTPUT){
        void *Vpt;
        Vpt = (void *)(Convert_dbl2flt(output->V[NVAR+ns],1.0, 0))[0][0];
        dataset = H5Dcreate(group, output->var_name[NVAR+ns], H5T_NATIVE_FLOAT,
                            dataspace, H5P_DEFAULT);

        #if MPI_POSIX == NO
        err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                       dataspace, plist_id_mpiio, Vpt);
        #else
        err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                       dataspace, H5P_DEFAULT, Vpt);
        #endif
      }

      H5Dclose(dataset);
      H5Sclose(memspace);
      H5Sclose(dataspace);
    }

    #if MPI_POSIX == NO
    H5Pclose(plist_id_mpiio);
    #endif
    H5Gclose(group);
  } /* Close group "stag_vars" */
#endif  /* STAGGERED_MHD */

  H5Gclose(timestep);

  group = H5Gcreate(file_identifier, "cell_coords", 0); /* Create group "cell_coords" (centered mesh) */ 
   
  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    dimens[nd] = grid->np_int_glob[nr];
  }
  dataspace = H5Screate_simple(rank, dimens, NULL);

#ifdef PARALLEL
  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    start[nd]  = grid->beg[nr] - grid->nghost[nr];
    stride[nd] = 1;
    count[nd]  = grid->np_int[nr];
  }
  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            start, stride, count, NULL);
#endif

  for (nd = 0; nd < DIMENSIONS; nd++){
    nr = DIMENSIONS-nd-1;
    dimens[nd] = grid->np_int[nr];
  }
  memspace = H5Screate_simple(rank,dimens,NULL);

/* ------------------------------------
       write cell centered mesh
   ------------------------------------ */

#if MPI_POSIX == NO
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
#endif

  for (nc = 0; nc < 3; nc++) {
    dataset = H5Dcreate(group, cname[nc], H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);

    #if MPI_POSIX == NO
    err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                   dataspace, plist_id_mpiio, cell_coords[nc][0][0]);
    #else
    err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                   dataspace, H5P_DEFAULT, cell_coords[nc][0][0]);
    #endif

    H5Dclose(dataset);
  }

#if MPI_POSIX == NO
  H5Pclose(plist_id_mpiio);
#endif
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Gclose(group); /* Close group "cell_coords" */

  group = H5Gcreate(file_identifier, "node_coords", 0); /* Create group "node_coords" (node mesh) */

  for (nd = 0; nd < DIMENSIONS; nd++){
    nr = DIMENSIONS-nd-1;
    dimens[nd] = grid->np_int_glob[nr] + 1;
  }

  dataspace = H5Screate_simple(rank, dimens, NULL);

#ifdef PARALLEL
  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    start[nd]  = grid->beg[nr] - grid->nghost[nr];
    stride[nd] = 1;
    count[nd]  = grid->np_int[nr];
    if (grid->rbound[nr] != 0) count[nd] += 1;
  }

  err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                            start, stride, count, NULL);
#endif

  for (nd = 0; nd < DIMENSIONS; nd++) {
    nr = DIMENSIONS-nd-1;
    dimens[nd] = grid->np_int[nr];
    if (grid->rbound[nr] != 0) dimens[nd] += 1;
  }
  memspace = H5Screate_simple(rank,dimens,NULL);

/* ------------------------------------
          write node centered mesh
   ------------------------------------ */

#if MPI_POSIX == NO
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
#endif

  for (nc = 0; nc < 3; nc++) {
    dataset = H5Dcreate(group, cname[nc], H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);

    #if MPI_POSIX == NO
    err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                   dataspace, plist_id_mpiio, node_coords[nc][0][0]);
    #else
    err = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace,
                   dataspace, H5P_DEFAULT, node_coords[nc][0][0]);
    #endif
    H5Dclose(dataset);
  }

#if MPI_POSIX == NO
  H5Pclose(plist_id_mpiio);
#endif
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Gclose(group); /* Close group "node_coords" */
  H5Fclose(file_identifier);

/* Create XDMF file to read HDF5 output (Visit or Paraview) */

  if (prank == 0) {

    if (output->type == DBL_H5_OUTPUT){
      sprintf(xmfext,"dbl.xmf");
      nprec = 8;
    } else if (output->type == FLT_H5_OUTPUT){
      sprintf(xmfext,"flt.xmf");
      nprec = 4;
    }
 
    sprintf (filenamexmf, "%s/data.%04d.%s", output->dir, output->nfile, xmfext);
   
    fxmf = fopen(filenamexmf, "w");
    fprintf(fxmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(fxmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(fxmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(fxmf, " <Domain>\n");
    fprintf(fxmf, "   <Grid Name=\"node_mesh\" GridType=\"Uniform\">\n");
    fprintf(fxmf, "    <Time Value=\"%12.6e\"/>\n",g_time);
    #if DIMENSIONS == 2
    fprintf(fxmf,"     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n",
            grid->np_int_glob[JDIR]+1, grid->np_int_glob[IDIR]+1);
    fprintf(fxmf, "     <Geometry GeometryType=\"X_Y\">\n");
    for (nd = 0; nd < DIMENSIONS; nd++) {
      fprintf(fxmf,"       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
             grid->np_int_glob[JDIR]+1, grid->np_int_glob[IDIR]+1 );
      fprintf(fxmf, "        %s:/node_coords/%s\n",filename,cname[nd]);
      fprintf(fxmf, "       </DataItem>\n");
    }
    #elif DIMENSIONS == 3
    fprintf(fxmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n",
            grid->np_int_glob[KDIR]+1,
            grid->np_int_glob[JDIR]+1,
            grid->np_int_glob[IDIR]+1);
    fprintf(fxmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    for (nd = 0; nd < DIMENSIONS; nd++) {
      fprintf(fxmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
              grid->np_int_glob[KDIR]+1,
              grid->np_int_glob[JDIR]+1,
              grid->np_int_glob[IDIR]+1 );
      fprintf(fxmf, "        %s:/node_coords/%s\n",filename,cname[nd]);
      fprintf(fxmf, "       </DataItem>\n");
    }
    #endif
    fprintf(fxmf, "     </Geometry>\n");
    for (nv = 0; nv < output->nvar; nv++) { /* Write cell-centered variables */
      if (!output->dump_var[nv] || output->stag_var[nv] != -1) continue; /* -- skip variable if excluded from output or if it is staggered */
      fprintf(fxmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
              output->var_name[nv]);
      #if DIMENSIONS == 2
      fprintf(fxmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
              grid->np_int_glob[JDIR],
              grid->np_int_glob[IDIR], nprec);
      #elif DIMENSIONS == 3
      fprintf(fxmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
              grid->np_int_glob[KDIR],
              grid->np_int_glob[JDIR],
              grid->np_int_glob[IDIR], nprec);
      #endif
      fprintf(fxmf, "        %s:%s/vars/%s\n",filename,tstepname,output->var_name[nv]);
      fprintf(fxmf, "       </DataItem>\n");
      fprintf(fxmf, "     </Attribute>\n");
    }
    for (nd = 0; nd < DIMENSIONS; nd++) { /* Write cell center coordinates (as cell-centerd variables) */
      fprintf(fxmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
              cname[nd]);
      #if DIMENSIONS == 2
      fprintf(fxmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
              grid->np_int_glob[JDIR], grid->np_int_glob[IDIR]);
      #elif DIMENSIONS == 3
      fprintf(fxmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
              grid->np_int_glob[KDIR],
              grid->np_int_glob[JDIR],
              grid->np_int_glob[IDIR]);
      #endif
      fprintf(fxmf, "        %s:/cell_coords/%s\n",filename,cname[nd]);
      fprintf(fxmf, "       </DataItem>\n");
      fprintf(fxmf, "     </Attribute>\n");
    } 
    fprintf(fxmf, "   </Grid>\n");
    fprintf(fxmf, " </Domain>\n");
    fprintf(fxmf, "</Xdmf>\n");
    fclose(fxmf);
  } /* if (prank == 0) */

/* XDMF file */

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  if (prank == 0){
    time(&tend);
    printLog (" [%5.2f sec]",difftime(tend,tbeg));
  }
#endif
}

/* ********************************************************************* */
void ReadHDF5 (Output *output, Grid *grid)
/*!
 *  Read data from disk using hdf5  format (double precision).
 *
 * \param [in] output the output structure associated with HDF5 format
 * \param [in] grid   a pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *********************************************************************** */
{
  hid_t dataspace, memspace, dataset;
  hid_t file_identifier, group, timestep;
  hid_t file_access = 0;
#if MPI_POSIX == NO
  hid_t plist_id_mpiio = 0; /* for collective MPI I/O */
#endif
  hid_t err;

  hsize_t dimens[DIMENSIONS];
  hsize_t start[DIMENSIONS];
  hsize_t stride[DIMENSIONS];
  hsize_t count[DIMENSIONS];

  char filename[128], tstepname[32];
  int ierr, rank, nd, nv, ns, nr;

/* --------------------------------------------------------------
     Since data is read in reverse order (Z-Y-X) it is
     convenient to define pointers to grid in reverse order
   -------------------------------------------------------------- */

  printLog ("> restarting from file #%d (dbl.h5)\n",output->nfile);
  sprintf (filename, "%s/data.%04d.dbl.h5", output->dir, output->nfile);

  rank = DIMENSIONS;

#ifdef PARALLEL
  file_access = H5Pcreate (H5P_FILE_ACCESS);
  #if MPI_POSIX == YES
  H5Pset_fapl_mpiposix(file_access, MPI_COMM_WORLD, 1);
  #else
  H5Pset_fapl_mpio(file_access,  MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  file_identifier = H5Fopen(filename, H5F_ACC_RDONLY, file_access);
  H5Pclose(file_access);
#else
  file_identifier = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
#endif

  if (file_identifier < 0){
    printLog ("! HDF5_READ: file %s does not exist\n");
    QUIT_PLUTO(1);
  }

  sprintf (tstepname, "Timestep_%d", output->nfile);
  timestep = H5Gopen(file_identifier, tstepname);
  group = H5Gopen(timestep, "vars");

#if MPI_POSIX == NO
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
#endif

  for (nv = 0; nv < output->nvar; nv++) {

  /* -- skip variable if excluded from output or if it is staggered -- */

    if (!output->dump_var[nv] || output->stag_var[nv] != -1) continue;

    dataset   = H5Dopen(group, output->var_name[nv]);
    dataspace = H5Dget_space(dataset);

    #ifdef PARALLEL
    for (nd = 0; nd < DIMENSIONS; nd++) {
      nr = DIMENSIONS-nd-1;
      start[nd]  = grid->beg[nr] - grid->nghost[nr];
      stride[nd] = 1;
      count[nd]  = grid->np_int[nr];
    }
    err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, NULL);  
    #endif

    for (nd = 0; nd < DIMENSIONS; nd++) {
      nr = DIMENSIONS-nd-1;
      dimens[nd] = grid->np_tot[nr];
    } 

    memspace = H5Screate_simple(rank,dimens,NULL);

    for (nd = 0; nd < DIMENSIONS; nd++){
      nr = DIMENSIONS-nd-1;
      start[nd]  = grid->nghost[nr];
      stride[nd] = 1;
      count[nd]  = grid->np_int[nr];
    }

    err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET, start, stride, count, NULL);

    #if MPI_POSIX == NO
    err = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                  plist_id_mpiio, output->V[nv][0][0]);
    #else
    err = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                  H5P_DEFAULT, output->V[nv][0][0]);
    #endif

    H5Dclose(dataset);
    H5Sclose(memspace);
    H5Sclose(dataspace);
  }

#if MPI_POSIX == NO
  H5Pclose(plist_id_mpiio);
#endif
  H5Gclose(group);

#ifdef STAGGERED_MHD
  group = H5Gopen(timestep, "stag_vars");

  #if MPI_POSIX == NO
  plist_id_mpiio = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id_mpiio,H5FD_MPIO_COLLECTIVE);
  #endif

  for (ns = 0; ns < DIMENSIONS; ns++) {
    dataset   = H5Dopen(group, output->var_name[NVAR+ns]);
    dataspace = H5Dget_space(dataset);

    #ifdef PARALLEL
    for (nd = 0; nd < DIMENSIONS; nd++) {
      nr = DIMENSIONS-nd-1;
      start[nd]  = grid->beg[nr] - grid->nghost[nr];
      stride[nd] = 1;
      count[nd]  = grid->np_int[nr];
      if (ns == DIMENSIONS-1-nd){
        if (grid->lbound[ns] != 0) count[nd] += 1;
        else                       start[nd] += 1;
      }
    }
    err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                              start, stride, count, NULL);
    #endif
 
    for (nd = 0; nd < DIMENSIONS; nd++){
      nr = DIMENSIONS-nd-1;
      dimens[nd] = grid->np_tot[nr] + (ns == nr);
    }

    memspace = H5Screate_simple(rank,dimens,NULL);

    for (nd = 0; nd < DIMENSIONS; nd++){
      nr = DIMENSIONS-nd-1;
      start[nd]  = grid->nghost[nr];
      stride[nd] = 1;
      count[nd]  = grid->np_int[nr];
      if (ns == nr && grid->lbound[ns] != 0) {
        start[nd] -= 1;
        count[nd] += 1;
      }
    }

    err = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,
                              start, stride, count, NULL); 
    #if MPI_POSIX == NO
    err = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace,
                  dataspace, plist_id_mpiio, output->V[NVAR+ns][0][0]);
    #else
    err = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace,
                  dataspace, H5P_DEFAULT, output->V[NVAR+ns][0][0]);
    #endif
    H5Dclose(dataset);
    H5Sclose(memspace); 
    H5Sclose(dataspace);
  }

#if MPI_POSIX == NO
  H5Pclose(plist_id_mpiio);
#endif
  H5Gclose(group);
#endif

  H5Gclose(timestep);
  H5Fclose(file_identifier);
}
