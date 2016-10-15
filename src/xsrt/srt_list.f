*+SRT_LIST	List all parameters
	SUBROUTINE SRT_LIST(IU,ISTAT)
	IMPLICIT NONE
	INTEGER IU,ISTAT
*IU	input	output channel number
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Dec-20
	INCLUDE 'SRT_COM'
	INTEGER J,I,K,NP,NDEF,NSQ,N1,N2,JJ
C
1000	FORMAT(1X,I5,'   type ',I3,'   def ',I2,'   sub ',I3,'   squa ',I3,
     +'   hit ',I5,'   miss ',I5)
1001	FORMAT(9X,'axis - surface normal      ',3(1X,G12.5),/,
     +        9X,'reference axis - tangent   ',3(1X,G12.5),/,
     +        9X,'vertex - reference position',3(1X,G12.5))
1002	FORMAT(1X,I5,'   submatrices ',I3,'   nx ',I5,'   ny ',I5)
1003	FORMAT(1X,I5,'   type ',I3)
1004	FORMAT(11X,5(1X,G12.5))
1005	FORMAT(
     +  1X,'Maximum number of parameters         ',I5,5X,I5,' set',/,
     +  1X,'Maximum number of deformations       ',I5,5X,I5,' set',/,
     +  1X,'Maximum number of deformation values ',I5,5X,I5,' set',/,
     +  1X,'Maximum number of surface qualities  ',I5,5X,I5,' set')
1006	FORMAT(1X,'   type ',I3,/,
     +        4X,'direction (for infinite source)',3(1X,G12.5),/,
     +        4X,'position (for finite source)   ',3(1X,G12.5),/,
     +        4X,'aperture area per ray (normal) ',1X,G12.5,/,
     +        4X,'normal to aperture             ',3(1X,G12.5),/,
     +        4X,'reference axis across aperture ',3(1X,G12.5),/,
     +        4X,'reference position in aperture ',3(1X,G12.5),/,
     +        4X,'limits ',4(1X,G12.5),/,
     +        4X,'deformation ',1X,I2,/,
     +        4X,'number of rays ',1X,I10)
1007	FORMAT(1X,'   surface ',I3)
1008	FORMAT(1X,'       XSAM(',I2,')')
1009	FORMAT(1X,'       YSAM(',I2,')')
1010	FORMAT(1X,'       ZDEF(',I2,')')
1011	FORMAT(1X,'       DZDX(',I2,')')
1012	FORMAT(1X,'       DZDY(',I2,')')
	IF(ISTAT.NE.0) RETURN
C
	WRITE(IU,*) 'Sequence of surface elements'
	DO J=1,NSUR
		WRITE(IU,1000) J,ISURS(J),ISDF(1,J),ISDF(2,J),ISTY(J),
     +		IHIT(J),IMISS(J)
		WRITE(IU,1001) (PAR(I),I=IPAR(J),IPAR(J)+8)
		DO K=IPAR(J)+9,IPAR(J)+MPAR(J)-1,5
			WRITE(IU,1004) (PAR(I),I=K,MIN(K+4,IPAR(J)+MPAR(J)-1))
		ENDDO
	ENDDO
	WRITE(IU,*) 'Deformations'
	NDEF=0
	DO J=1,MAXDF
		IF(IDFM(1,J).NE.0) THEN
			NDEF=NDEF+1
			WRITE(IU,1002) J,IDFM(1,J),IDFM(2,J),IDFM(3,J)
      			DO JJ=1,IDFM(1,J)
     				WRITE(IU,1008) JJ
				N1=IDFP(1,J)+(JJ-1)*IDFM(2,J)
      				N2=N1+IDFM(2,J)-1
				DO K=N1,N2,5
 				 WRITE(IU,1004) (DSAM(I),I=K,MIN(K+4,N2))
				ENDDO
     				WRITE(IU,1009) JJ
				N1=IDFP(2,J)+(JJ-1)*IDFM(3,J)
      				N2=N1+IDFM(3,J)-1
				DO K=N1,N2,5
 				 WRITE(IU,1004) (DSAM(I),I=K,MIN(K+4,N2))
				ENDDO
     				WRITE(IU,1010) JJ
				N1=IDFP(3,J)+(JJ-1)*IDFM(2,J)*IDFM(3,J)
      				N2=N1+IDFM(2,J)*IDFM(3,J)-1
				DO K=N1,N2,5
 				 WRITE(IU,1004) (DSAM(I),I=K,MIN(K+4,N2))
				ENDDO
     				WRITE(IU,1011) JJ
				N1=IDFP(4,J)+(JJ-1)*IDFM(2,J)*IDFM(3,J)
      				N2=N1+IDFM(2,J)*IDFM(3,J)-1
				DO K=N1,N2,5
 				 WRITE(IU,1004) (DSAM(I),I=K,MIN(K+4,N2))
				ENDDO
     				WRITE(IU,1012) JJ
				N1=IDFP(5,J)+(JJ-1)*IDFM(2,J)*IDFM(3,J)
      				N2=N1+IDFM(2,J)*IDFM(3,J)-1
				DO K=N1,N2,5
 				 WRITE(IU,1004) (DSAM(I),I=K,MIN(K+4,N2))
				ENDDO
			ENDDO
		ENDIF
	ENDDO
	WRITE(IU,*) 'Surface quality'
	NSQ=0
	DO J=1,MAXST
		IF(ISQP(1,J).NE.0) THEN
			NSQ=NSQ+1
			WRITE(IU,1003) J,ISQP(1,J)
			DO K=ISQP(2,J),ISQP(2,J)+ISQP(3,J)-1,5
				WRITE(IU,1004) (PAR(I),I=K,MIN(K+4,
     +				ISQP(2,J)+ISQP(3,J)-1))
			ENDDO
		ENDIF
	ENDDO
	WRITE(IU,*) 'Source'
	IF(ISRC(1).GT.0) THEN
	  WRITE(IU,1006) ISRC(1),(PAR(I),I=ISRC(2),ISRC(2)+19),
     +		ISRC(3),NRAYS
	ENDIF
	WRITE(IU,*) 'Detector'
	IF(IDET.GT.0) THEN
	  WRITE(IU,1007) IDET
	ENDIF
	WRITE(IU,1005) MAXPAR,NPAR,MAXDF,NDEF,MAXDSAM,IDSAM,MAXST,NSQ
	END
