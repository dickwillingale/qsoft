/*
  Fortran interface for C routines on ismatten.c
  Dick Willingale 1994-Nov-2
*/
double atten_(wav,hcol,heicol,heiicol)
double *wav,*hcol,*heicol,*heiicol;
{
    extern double atten();
    double trans;
    trans = atten(*wav,*hcol,*heicol,*heiicol);
    return(trans);
}
double tauh_(wav,hcol,zee)
double *wav,*hcol,*zee;
{
    extern double tauh();
    double tau;
    tau=tauh(*wav,*hcol,*zee);
    return(tau);
}
double tauhe_(wav,hecol)
    double *wav,*hecol;
{
    extern double tauhe();
    double tau;
    tau=tauhe(*wav,*hecol);
    return(tau);
}
