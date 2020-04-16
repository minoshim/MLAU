void bkup_save(int n, int cnt, double tim, double dt, double trec, int mpi_rank)
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
  fwrite(ro,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_mx_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(mx,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_my_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(my,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_mz_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(mz,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_bx_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(bx,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_by_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(by,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_bz_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(bz,sizeof(double),nx*ny,outfil);
  fclose(outfil);
  sprintf(filname,"%s/bkup_en_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  fwrite(en,sizeof(double),nx*ny,outfil);
  fclose(outfil);
}
