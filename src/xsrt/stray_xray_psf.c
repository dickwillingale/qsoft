#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.141592653589793115997963468544
#define MMAX 1100
#define EMAX 10
#define AMAX 20

double stray_xray_sb(double ekev,double phi,int na,int ne,
	double rad[na],double kev[ne],double sb[EMAX][AMAX]) {
/* Surface brightness of stray x-ray PSF at r=0
   ekev		input	photon energy keV
   phi		input	source off-axis angle, radians
   na		input	number of angle samples
   ne		input	numner of energy samples
   rad		input	off-axis angles, arcmins
   kev		input	photon energies, keV
   sb		input	surface brightness at r=0 (source position) cm2/arcmins2
   return surface brightness at r=0, cm2/arcmins2
   linear interpolation from array sb
   The coefficients in array sb estimated by raytracing with an appropriate
   Wolter I aperture geometry and reflecting surface coating.
   e.g.  Athena SPO modules Ir + B4C, Swift XRT Au

   Dick Willingale 2016-Apr-17
*/
	double am,sbl,sbh,rat;
	int ia,ie,i;
	am=phi*180.0*60.0/PI;
/* Find indices in array */
	ia=0;
	for(i=1;i<(na-1);i++) {
		if(am>rad[i]) ia=i;
	}
	ie=0;
	for(i=1;i<(ne-1);i++) {
		if(ekev>kev[i]) ie=i;
	}
/* Linear interpolation in angle at 2 energies */
	rat=(am-rad[ia])/(rad[ia+1]-rad[ia]);
	sbl=sb[ie][ia]+(sb[ie][ia+1]-sb[ie][ia])*rat;
	sbh=sb[ie+1][ia]+(sb[ie+1][ia+1]-sb[ie+1][ia])*rat;
/* Linear interpolation between 2 energies */
	rat=(ekev-kev[ie])/(kev[ie+1]-kev[ie]);
	return(sbl+(sbh-sbl)*rat);
}

double stray_xray_psf(double ekev,double x,double y,double xs,double ys,
       	int nm, double tg1[nm],double tg2[nm],
       	double th1[nm], double th2[nm],
	int na,int ne,double rad[na],double kev[ne],double sb[EMAX][AMAX]) {
/* Stray X-ray PSF surface brightness from Wolter I aperture array at specified
   position in FOV
   ekev	input	photon energy keV
   x	input	azimuthal position in FOV, radians
   y	input	elevation position in FOV, radians
   xs	input	azimuth of source outside FOV, radians
   ys	input	elevation of source outside FOV, radians
   nm	input	number of Wolter I aperture sectors
   tg1	input	minimum on-axis grazing angle for sector radians
   tg2	input	maximum on-axis grazing angle for sector, radians
   th1	input	minimum angular position of sector aperture, radians
   th2	input	maximum angular position of sector aperture, radians
   na	input	number of radial samples for normalisation
   ne	input	number of energy samples for normalisation
   rad	input	off-axis angles, arcmins
   kev	input	photon energies, keV
   sb	input	surface brightness at r=0 (source position) cm2/arcmins2
   return surface brightness of stray PSF at x,y, cm2/arcmins2

   Dick Willingale 2016-Apr-17
*/
	double phi,alpha,r1,r2,t1,t2,th,r,xx,yy,pc,rmax,b;
	int im;
/* phi is off-axis angle of the source, alpha is the angular position of the
   source wrt the x axis. i.e. phi, alpha give the position of the source in
   polar coordinates */
	phi=sqrt(xs*xs+ys*ys);
	alpha=atan2(ys,xs);
/* rotate the coordinate system so that source lies at origin on the x-axis */
	xx= phi-x*cos(alpha)-y*sin(alpha);
	yy= -x*sin(alpha)+y*cos(alpha);
/* calculate the polar coordiates about the source position */
	th=atan2(yy,xx);
	r=sqrt(xx*xx+yy*yy);
/* calculate the maximum radius from the source and generate a linear
   profile for the surface brightness */
	pc=phi*cos(th);
	rmax=2.0*(pc-phi*0.5);
	b=(rmax-r)/rmax;
/* Loop through aperture sectors to see if the stray patch from any sector
   overlaps with the position in the FOV */
	if(b>0) {
	    for(im=0;im<nm;im++) {
		t2=alpha-th1[im]+PI;
		if(t2>PI) {
			t2=t2-PI*2.0;
		}
		t1=alpha-th2[im]+PI;
		if(t1>PI) {
			t1=t1-PI*2.0;
		}
		r1=2.0*(pc-tg1[im]);
		r2=2.0*(pc-tg2[im]);
		if(r<r1 && r>r2 && th>t1 && th<t2) {
/* overlap found - set surface brightness and return */
			b=b*stray_xray_sb(ekev,phi,na,ne,rad,kev,sb);
			return(b);
		}
	    }
	}
	b=0.0;
	return(b);
}
void stray_xray_image(double *ekev,
		int *nx, double x[*nx], int *ny, double y[*ny],
		double *xs, double *ys,double arr[*ny][*nx])
{
/* Generate image of stray X-rays from Wolter I aperture sector array
   ekev	input	photon energy
   nx	input	number of x grid points (columns in image)
   x	input	azimuthal positions of x grid centre of pixels, radians
   ny	input	number of y grid points (rows in image)
   y	input	elevation positions of y grid centre of pixels, radians
   xs	input	azimuth of source, radians
   ys	input	elevation of source, radians
   arr	output	image of PSF cm^2/arcmins^2

   Dick Willingale 2016-Apr-17
*/
	double tg1[MMAX],tg2[MMAX],th1[MMAX],th2[MMAX];
	double rad[AMAX],kev[EMAX],sb[EMAX][AMAX];
	int ix,iy,nm,na,ne,im;
	char header[80];
	FILE *fh;
/* read in Wolter I aperture sector data from file */
	fh = fopen("wolter_1_aperture_sectors.dat","r");
	if(fh==NULL) {
		printf("Cannot open file wolter_1_aperture_sectors.dat\n");
		exit(EXIT_FAILURE);
	}
	fgets(header,80,fh);
	nm=0;
	while(!feof(fh) && nm<MMAX) {
		if( fscanf(fh,"%d %lf %lf %lf %lf", &im,
		&tg1[nm],&tg2[nm],&th1[nm],&th2[nm])>0) nm++;
	}
	fclose(fh);
	fh = fopen("wolter_1_stray_psf.dat","r");
	if(fh==NULL) {
		printf("Cannot open file wolter_1_stray_psf.dat\n");
		exit(EXIT_FAILURE);
	}
	ne=4;
	fscanf(fh,"%lf %lf %lf %lf", &kev[0],&kev[1],&kev[2],&kev[3]);
	na=0;
	while(!feof(fh) && na<AMAX) {
		if( fscanf(fh,"%lf %lf %lf %lf %lf",
		 &rad[na],&sb[0][na],&sb[1][na],
		 &sb[2][na],&sb[3][na])>0) na++;
	}
	fclose(fh);
/* loop for all pixels in image */
	for(ix=0;ix<*nx;ix++) {
		for(iy=0;iy<*ny;iy++) {
			arr[iy][ix]=stray_xray_psf(*ekev,x[ix],y[iy],
			*xs,*ys,nm,tg1,tg2,th1,th2,na,ne,rad,kev,sb);
		}
	}
}
