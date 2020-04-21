void dataio(int n, int cnt, double tim)
// Output simulation data
{
  char filname[100];
  FILE *outfil;

  if (n == 0){
    sprintf(filname,"%s/params.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%f\n",gam);
    fclose(outfil);

    sprintf(filname,"%s/offsets.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%d %d %d\n",XOFF,YOFF,ZOFF);
    fclose(outfil);
    
    sprintf(filname,"%s/t.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%f\n",tim);
    fclose(outfil);
  } else{
    sprintf(filname,"%s/t.dat",fildir);
    outfil=fopen(filname,"a");
    fprintf(outfil,"%f\n",tim);
    fclose(outfil);
  }

  sprintf(filname,"%s/outdat_%05d.dat",fildir,cnt);  
  outfil=fopen(filname,"wb");
  fwrite(&ro[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&mx[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&my[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&mz[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&en[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&bx[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&by[0],sizeof(double),nx*ny*nz,outfil);
  fwrite(&bz[0],sizeof(double),nx*ny*nz,outfil);
  fclose(outfil);

}
