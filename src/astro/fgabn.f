*+FGABN get element abundance
	subroutine fgabn(celt,abu)
        implicit none
	character*2 celt
        real abu
*- Dick Willingale 2014-Apr-08
	include 'SPX_COM'
	integer i
	do i=1,nels
	     if(element(i) .eq. celt) then
		     iels=i
                     abu=abund(iabu,iels)
		     if(i.gt.2) then
			     abu=abu*metalicity
		     endif
		     return
	     endif
	enddo
	write(*,*) 'fgabn error: element ',celt,' not found'
	end
