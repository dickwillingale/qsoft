*+SRT_SETF	Set surface form and limits parameters
	SUBROUTINE SRT_SETF(NS,IT,NP,P,IDEF,IQ,IH,IM,ISTAT)
	IMPLICIT NONE
	INTEGER NS,IT,NP,IDEF(2),IQ,IH,IM,ISTAT
	DOUBLE PRECISION P(NP)
*NS	input	surface number (0 for new entry)
*IT	input	surface type
*NP	input	number of parameters
*P	input	array of parameters
*IDEF	input	deformation
*IQ	input	surface quality
*IH	input	hit index (-ve for next in sequence)
*IM	input	miss index (-ve for next in sequence)
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Dec-6
	INCLUDE 'SRT_COM'
	INTEGER NSI,NPI,J
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NS.EQ.0) THEN
C New surface
		IF(NSUR.EQ.MAXSUR) THEN
			WRITE(*,*) 'SRT_SETF error - surface list full'
			WRITE(*,*) 'MAXSUR ',MAXSUR
			WRITE(*,*) 'MAXPAR ',MAXPAR
			WRITE(*,*) 'requires NPAR ',NPAR+NP,' NSUR ',NSUR+1
			ISTAT=1
			RETURN
		ENDIF
		NSUR=NSUR+1
		NSI=NSUR
C Check parameter space
		IF(NPAR+NP.GT.MAXPAR) THEN
			WRITE(*,*) 'SRT_SETF error - parameter space full'
			WRITE(*,*) 'MAXSUR ',MAXSUR
			WRITE(*,*) 'MAXPAR ',MAXPAR
			WRITE(*,*) 'requires NPAR ',NPAR+NP,' NSUR ',NSUR+1
			ISTAT=1
			RETURN
		ENDIF
C Allocate parameter index
		NPI=NPAR+1
		IPAR(NSI)=NPI
		NPAR=NPAR+NP
	ELSE
		IF(NS.GT.NSUR) THEN
			WRITE(*,*) 
     +			'SRT_SETF error - attempt to modify null surface'
                        WRITE(*,*) 'NS NSUR',NS,NSUR
			ISTAT=1
			RETURN
		ENDIF
C Get surface and parameter index
		NSI=NS
		NPI=IPAR(NSI)
	ENDIF
C Proceed with setting
	DO J=1,NP
		PAR(NPI+J-1)=P(J)	
	ENDDO
	MPAR(NSI)=NP
	ISURS(NSI)=IT
	ISDF(1,NSI)=IDEF(1)
	ISDF(2,NSI)=IDEF(2)
	ISTY(NSI)=IQ
	IF(IH.EQ.-1) THEN
		IHIT(NSI)=NSI+1
	ELSEIF(IH.EQ.-2) THEN
		IHIT(NSI)=NSI
	ELSE	
		IHIT(NSI)=IH
	ENDIF
	IF(IM.EQ.-1) THEN
		IMISS(NSI)=NSI+1
	ELSEIF(IM.EQ.-2) THEN
		IMISS(NSI)=NSI
	ELSE
		IMISS(NSI)=IM
	ENDIF
C If a detecting surface set detector surface index
	IF(IQ.LT.0) THEN
		IDET=NSI
	ENDIF
	END
