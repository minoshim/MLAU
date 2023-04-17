inline void conv_d2f(double vali[], float valo[], unsigned long nn)
// convert double to single precision
{
  unsigned long i;
  for (i=0;i<nn;i++) valo[i]=(float)vali[i];
}

void dataio(int n, int cnt, double tim, int mpi_rank)
{
  int i;
  int nall=nx*ny;
  FILE *outfil;
  char filname[100];
  float *val;
  val=new float[nall];

  if (n == 0){
    if (mpi_rank == 0){
      sprintf(filname,"%s/params.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%f %f %f\n",gam,beta,lambda);
      fclose(outfil);
      
      sprintf(filname,"%s/offsets.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%d %d\n",xoff,yoff);
      fclose(outfil);

      sprintf(filname,"%s/t.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%f\n",tim);
      fclose(outfil);
    }
    if ((mpi_rank/mpi_numx) == 0){
      sprintf(filname,"%s/x_%05d.dat",fildir,(mpi_rank%mpi_numx));
      outfil=fopen(filname,"w");
      for (i=0;i<nx;i++) fprintf(outfil,"%.12f\n",x[i]);
      fclose(outfil);
    }
    if ((mpi_rank%mpi_numx) == 0){
      sprintf(filname,"%s/y_%05d.dat",fildir,(mpi_rank/mpi_numx));
      outfil=fopen(filname,"w");
      for (i=0;i<ny;i++) fprintf(outfil,"%.12f\n",y[i]);
      fclose(outfil);
    }
  } else{
    if (mpi_rank == 0){
      sprintf(filname,"%s/t.dat",fildir);
      outfil=fopen(filname,"a");
      fprintf(outfil,"%f\n",tim);
      fclose(outfil);
    }
  }

  sprintf(filname,"%s/outdat_%05d_%05d.dat",fildir,cnt,mpi_rank);
  outfil=fopen(filname,"wb");
  conv_d2f(ro,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(mx,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(my,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(mz,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(en,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(bx,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(by,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  conv_d2f(bz,val,nall);
  fwrite(val,sizeof(float),nall,outfil);
  fclose(outfil);

  delete[] val;
}
