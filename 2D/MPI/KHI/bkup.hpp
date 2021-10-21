void bkup_load(int *n, int *cnt, double *tim, double *dt, double *trec, int nall, int mpi_rank);
void bkup_save(int n, int cnt, double tim, double dt, double trec, int nall, int mpi_rank);

void bkup_load(int *n, int *cnt, double *tim, double *dt, double *trec, int nall, int mpi_rank)
// Load bkup files
{
  FILE *infil;
  int vali[2];
  double vald[3];
  char filname[100];
  int itmp;
  size_t ttmp;
  
  sprintf(filname,"%s/bkup_stamp.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    itmp=fscanf(infil,"%d %d\n",&vali[0],&vali[1]);
    fclose(infil);
    *n=vali[0];
    *cnt=vali[1];

    sprintf(filname,"%s/bkup_time.dat",fildir);
    infil=fopen(filname,"rb");
    ttmp=fread(&vald,sizeof(double),3,infil);
    fclose(infil);
    *tim=vald[0];
    *dt=vald[1];
    *trec=vald[2];
    
    sprintf(filname,"%s/bkup_ro_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(ro,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_mx_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(mx,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_my_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(my,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_mz_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(mz,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_bx_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(bx,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_by_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(by,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_bz_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(bz,sizeof(double),nall,infil);
    fclose(infil);
    sprintf(filname,"%s/bkup_en_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    ttmp=fread(en,sizeof(double),nall,infil);
    fclose(infil);

    if (mpi_rank == 0) printf("Load backup files at %d steps (T = %f)\n",vali[0],vald[0]);
  }

}

void bkup_save(int n, int cnt, double tim, double dt, double trec, int nall, int mpi_rank)
// Save bkup files
{
  FILE *outfil;
  double vald[]={tim,dt,trec};
  char filname[100];

  if (mpi_rank == 0){
    sprintf(filname,"%s/bkup_stamp.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%d %d\n",n,cnt);
    fclose(outfil);

    sprintf(filname,"%s/bkup_time.dat",fildir);
    outfil=fopen(filname,"wb");
    fwrite(&vald,sizeof(double),3,outfil);
    fclose(outfil);
  }

  sprintf(filname,"%s/bkup_ro_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(ro,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_mx_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(mx,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_my_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(my,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_mz_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(mz,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_bx_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(bx,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_by_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(by,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_bz_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(bz,sizeof(double),nall,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_en_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(en,sizeof(double),nall,outfil);
  fclose(outfil);
}
