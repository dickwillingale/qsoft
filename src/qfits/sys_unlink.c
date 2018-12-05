/* SYS_UNLINK(FILENAME,ISTAT)		 */
/* CHARACTER*(*) FILENAME	input	 */
/* INTEGER ISTAT		in/out	 */
/* Fortran callable routine to delete a file */
#include <unistd.h>

void sys_unlink_(char *filename, int *istat, int lfile)

{
char* q=filename;
char name[256] ;
char* p=name;
int i;
for(i=0;i<lfile;i++) *p++=*q++;
*p=0;
*istat=unlink(name) ;
}
