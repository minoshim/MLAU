void init_grid()
// Define X and Y coordinates
{
  char filname[100];
  FILE *outfil;

  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (int i=0;i<nx;i++){
    // x[i]=(i-XOFF+0.5)*dx;
    x[i]=(i-XOFF)*dx;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (int j=0;j<ny;j++){
    // y[j]=(j-YOFF+0.5)*dy;
    y[j]=(j-YOFF)*dy;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}
