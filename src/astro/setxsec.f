*-SETXSEC	Set atomic cross-section data type for SPX
	subroutine setatomic(xsec)
        implicit none
        character*4 xsec
	real metal
*xsec      input     X-section - vern, bcmc obcm
*-Author Dick Willingale 2014-Apr-08
	include 'SPX_COM'
	integer i
C 
	if(xsec.eq.'vern'.or.xsec.eq.'bcmc'.or.xsec.eq.'obcm') then
	 	xsection=xsec
	else
		write(*,*) 'setxsec error: unknown xsec',xsec
	endif
	end
