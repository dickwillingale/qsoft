*+SRT_SAVE	Read rays from file into arrays
	SUBROUTINE SRT_SAVE(NP,RPOS,IQU,NR,RDIR,RQU,ND,DNML,DRXL
     +	,DPOS,RIRIS,
     +	IADJUST,PSAS,PPLANE,NDIFF,IM,DIM,Y,Z,WL,ND3,WK,
     +	DSHFT,YBAR,ZBAR,RMS,GY,GZ,
     +	XR,YR,ZR,XC,YC,ZC,XD,YD,ZD,AD,DL,IR,IB,ID,ISTAT)
	IMPLICIT NONE
	INTEGER NP,NR,IQU(NP),ND,IR(ND),IB(ND),ID(ND),ISTAT
	INTEGER NDIFF,ND3
	DOUBLE PRECISION RPOS(3,NP),RDIR(3,NR),RQU(3,NR)
	DOUBLE PRECISION XR(ND),YR(ND),ZR(ND),XC(ND),YC(ND),ZC(ND)
	DOUBLE PRECISION XD(ND),YD(ND),ZD(ND),AD(ND),DL(ND)
	DOUBLE PRECISION Y(ND),Z(ND),WL(ND)
	DOUBLE PRECISION DNML(3),DRXL(3),DPOS(3),DSHFT,YBAR,ZBAR,RMS,RIRIS
	DOUBLE PRECISION GY(ND),GZ(ND),PSAS,IM(NDIFF,NDIFF),WK(ND3)
	DOUBLE PRECISION DIM(NDIFF,NDIFF),PPLANE
	INTEGER IADJUST
*NP	input	maximum total number of positions on rays
*RPOS	output	positions on rays
*IQU	output	history of positions (-2 for start, 0 or -1 for end)
*NR	input	number of rays
*RDIR	output	direction of rays
*RQU	output	quality of rays
*ND	input	maximum number of detections
*DNML	input	normal to detector surface 
*DRXL	input	reference axis in surface of detector
*DPOS	input	origin of detector surface
*RIRIS	input	radius of iris for blur calculation
*IADJUSTinput	0 no adjustment, 1 x adjustment/shift using rms radius
*PSAS	input	pixel size arc seconds for diffraction (0.0 for no calculation)
*PPLANE	input	position of principal plane wrt detector
*NDIFF	input	size of arrays for diffraction calculation
*IM	input	work array for diffraction calculation
*DIM	output	diffraction limited image array (pixel size PSAS)
*Y	output	detected positions in detector frame
*Z	output	detected positions in the detector frame
*WL	output	work array
*ND3	input	size of work array
*WK	output	work array
*DSHFT 	output	detector axial shift
*YBAR	output	mean detector y position
*ZBAR	output	mean detector z position
*GY	output	gradient of ray wrt detector shift in x
*GZ	output	gradient of ray wrt detector shift in x
*RMS	output	root mean square radius of detected positions
*XR	output	last x ray position before detection
*YR	output	last y ray position before detection
*ZR	output	last z ray position before detection
*XC	output	ray direction x
*YC	output	ray direction y
*ZC	output	ray direction z
*XD	output	detected x positions
*YD	output	detected y positions
*ZD	output	detected z positions
*AD	output	areas detected
*DL	output	distance along ray to detector
*IR	output	number of hits surface index 1
*IB	output	number of hits surface index 2
*ID	output	number of hits surface index 3
*ISTAT	in/out	returned status
* Note the number of positions, rays and detections on the file
* is known in advance so checking is a precaution.
*-Author Dick Willingale 1996-Dec-10
	INCLUDE 'SRT_COM'
	INTEGER NPR,K,J,NC,I,KD,IY,IZ,IYY,IZZ,NELS(2),NN,IFAIL,JJ,KK,ND2
	DOUBLE PRECISION YAX(3),X,A,PP,DD,YPP,ZPP,PHY,PHZ,PSRAD,SI
	DOUBLE PRECISION DR(3),YB,ZB,Y2,Z2,YPB,ZPB,YYB,ZZB,YP2,ZP2
	DOUBLE PRECISION GX,RR,SS
	DOUBLE PRECISION DX,DY,DZ
C
	IF(ISTAT.NE.0) RETURN
C
	KD=0
	NC=1
      	REWIND ITRA
	DO K=1,NR
C Get ray from file
                READ(ITRA,IOSTAT=ISTAT) 
     +          NPR,((RPOS(I,J),I=1,3),J=NC,NPR+NC-1),(RDIR(I,K),I=1,3),
     +		(RQU(I,K),I=1,3),(IQU(J),J=NC,NPR+NC-1)
                IF(ISTAT.NE.0) THEN
                        WRITE(*,*) 'SRT_SAVE error - end of file'
			RETURN
                ENDIF
		NC=NC+NPR
		IF(NC.GT.NP+1) THEN
			WRITE(*,*) 'SRT_SAVE error - too many positions'
			ISTAT=1
			RETURN
		ENDIF
		IF(IQU(NC-1).EQ.-1) THEN
C Detected ray
			KD=KD+1
			IF(KD.GT.ND) THEN
				WRITE(*,*)
     +				'SRT_SAVE error - too many detections'
				ISTAT=1
				RETURN
			ENDIF
			XR(KD)=RPOS(1,NC-2)
			YR(KD)=RPOS(2,NC-2)
			ZR(KD)=RPOS(3,NC-2)
			XD(KD)=RPOS(1,NC-1)
			YD(KD)=RPOS(2,NC-1)
			ZD(KD)=RPOS(3,NC-1)
			XC(KD)=RDIR(1,K)
			YC(KD)=RDIR(2,K)
			ZC(KD)=RDIR(3,K)
			WL(KD)=RQU(1,K)
			AD(KD)=RQU(2,K)
C Count up history of ray
			IR(KD)=0
			IB(KD)=0
			ID(KD)=0
			DO I=NC-NPR+1,NC-2
				IF(IQU(I).EQ.1) THEN
					IR(KD)=IR(KD)+1
				ELSEIF(IQU(I).EQ.2) THEN
					IB(KD)=IB(KD)+1
				ELSEIF(IQU(I).EQ.3) THEN
					ID(KD)=ID(KD)+1
				ENDIF
			ENDDO
C Find length along ray
			DL(KD)=0.0
			DO I=NC-NPR+1,NC-1
				DX=RPOS(1,I)-RPOS(1,I-1)
				DY=RPOS(2,I)-RPOS(2,I-1)
				DZ=RPOS(3,I)-RPOS(3,I-1)
				DL(KD)=DL(KD)+SQRT(DX*DX+DY*DY+DZ*DZ)
			ENDDO
		ENDIF
	ENDDO
	IF(KD.GT.0) THEN
C Find third detector axis
		CALL SRT_VCRS(DNML,DRXL,YAX)
		CALL SRT_VNRM(YAX,ISTAT)
C Calculate area weighted mean and rms radius of detected positions
		YB=0.0
		ZB=0.0
		Y2=0.0
		Z2=0.0
		YPB=0.0
		ZPB=0.0
		YYB=0.0
		ZZB=0.0
		YP2=0.0
		ZP2=0.0
		SS=0.0
		DO J=1,KD
C Find position in detector coordinates
			DR(1)=XD(J)-DPOS(1)
			DR(2)=YD(J)-DPOS(2)
			DR(3)=ZD(J)-DPOS(3)
			CALL SRT_VDOT(DNML,DR,X)
			CALL SRT_VDOT(DRXL,DR,Y(J))
			CALL SRT_VDOT(YAX,DR,Z(J))
C Find gradients with respect to x-shift
			DR(1)=XC(J)
			DR(2)=YC(J)
			DR(3)=ZC(J)
			CALL SRT_VDOT(DR,DNML,GX)
			IF(GX.NE.0.0) THEN
				CALL SRT_VDOT(DR,DRXL,GY(J))
				CALL SRT_VDOT(DR,YAX,GZ(J))
				GY(J)=GY(J)/GX
				GZ(J)=GZ(J)/GX
			ELSE
				GY(J)=0.0
				GZ(J)=0.0
			ENDIF
C Only rays with IADJUST reflections within iris radius count
			IF(SQRT(Y(J)**2+Z(J)**2).LT.RIRIS.AND.
     +			IR(J).EQ.IADJUST) THEN
				A=AD(J)
				SS=SS+A
				YB=YB+Y(J)*A
				ZB=ZB+Z(J)*A
				Y2=Y2+A*Y(J)**2
				Z2=Z2+A*Z(J)**2
				YPB=YPB+GY(J)*A
				ZPB=YPB+GZ(J)*A
				YYB=YYB+Y(J)*GY(J)*A
				ZZB=ZZB+Z(J)*GZ(J)*A
				YP2=YP2+GY(J)*GY(J)*A
				ZP2=ZP2+GZ(J)*GZ(J)*A
			ENDIF
		ENDDO
		IF(SS.GT.0.0) THEN
			YB=YB/SS
			ZB=ZB/SS
			Y2=Y2/SS
			Z2=Z2/SS
			YPB=YPB/SS
			ZPB=ZPB/SS
			YYB=YYB/SS
			ZZB=ZZB/SS
			YP2=YP2/SS
			ZP2=ZP2/SS
			RR=YP2+ZP2-YPB**2-ZPB**2
		ELSE
			RR=0.0
		ENDIF
		IF(IADJUST.GT.0.AND.RR.NE.0.0) THEN
			IF(IADJUST.GT.0) THEN
C Shift detector along normal to position of minimum rms radius
				DSHFT=(YPB*YB+ZPB*ZB-YYB-ZZB)/RR
C			ELSEIF(IADJUST.EQ.2) THEN
C Shift detector along normal to position of minimum rms y
C				DSHFT=(YPB*YB-YYB)/(YP2-YPB**2)
C			ELSE
C Shift detector along normal to position of minimum rms z
C				DSHFT=(ZPB*ZB-ZZB)/(ZP2-ZPB**2)
			ENDIF
			DR(1)=DSHFT
			DO J=1,KD
				DR(2)=DSHFT*GY(J)
				DR(3)=DSHFT*GZ(J)
				Y(J)=Y(J)+DR(2)
				Z(J)=Z(J)+DR(3)
				XD(J)=XD(J)+DR(1)*DNML(1)+
     +				DR(2)*DRXL(1)+DR(3)*YAX(1)
				YD(J)=YD(J)+DR(1)*DNML(2)+
     +				DR(2)*DRXL(2)+DR(3)*YAX(2)
				ZD(J)=ZD(J)+DR(1)*DNML(3)+
     +				DR(2)*DRXL(3)+DR(3)*YAX(3)
			ENDDO
		ELSE
			DSHFT=0.0
		ENDIF
C Return mean and rms including effect of any shift
		YBAR=DSHFT*YPB+YB
		ZBAR=DSHFT*ZPB+ZB
		RMS=(YP2+ZP2-YPB**2-ZPB**2)*DSHFT**2+
     +		(YYB+ZZB-YPB*YB-ZPB*ZB)*DSHFT*2.D0+Y2+Z2-YB**2-ZB**2
		RMS=SQRT(ABS(RMS))
		IF(PSAS.GT.0.0) THEN
C Calculate diffraction limited image
C Find position of principal plane wrt shifted detector
			PP=PPLANE-DSHFT
C Calculate sample size in image plane
			PSRAD=PSAS*PI/3600/180
			SI=PSRAD*PP
C Calculate scaling factor for samples in principal plane
C Includes conversion from Angstroms to mm for wavelength
			DD=1.D-7/(DBLE(NDIFF)*PSRAD)
C Initialize real and imaginary part of principal plane arrays
			DO J=1,NDIFF
				DO K=1,NDIFF
					DIM(K,J)=0.0
					IM(K,J)=0.0
				ENDDO
			ENDDO
C Set maximum image sample index
			ND2=NDIFF/2
C Loop for all rays to bin up rays in principal plane
			DO J=1,KD
			    IF(WL(J).GT.0.0) THEN
C Calculate position of ray on principal plane
				YPP=Y(J)+PPLANE*GY(J)
				ZPP=Z(J)+PPLANE*GZ(J)
C Calculate sample position in the principal plane (note aliasing)
				IF(YPP.GT.0.0) THEN
				    IYY=MOD(INT(YPP/(WL(J)*DD)),NDIFF)+1
				ELSE
				    IYY=NDIFF-MOD(INT(-YPP/(WL(J)*DD)),NDIFF)
				ENDIF
				IF(ZPP.GT.0.0) THEN
				    IZZ=MOD(INT(ZPP/(WL(J)*DD)),NDIFF)+1
				ELSE
				    IZZ=NDIFF-MOD(INT(-ZPP/(WL(J)*DD)),NDIFF)
				ENDIF
C Calculate sample position in image plane
				IF(Y(J).GT.0.0) THEN
					IY=INT(Y(J)/SI)+1
				ELSE
					IY=INT(Y(J)/SI)
				ENDIF
				IF(Z(J).GT.0.0) THEN
					IZ=INT(Z(J)/SI)+1
				ELSE
					IZ=INT(Z(J)/SI)
				ENDIF
				IF(ABS(IY).LT.ND2.AND.ABS(IZ).LT.ND2) THEN
C Only include rays within the limits imposed by image sample range
C Calculate phase errors (note wavelength is in Angstoms, convert to mm)
					PHY=-2.D7*PI*IY*IYY*SI*DD/PP
					PHZ=-2.D7*PI*IZ*IZZ*SI*DD/PP
C Add complex amplitude of ray to principal plane
C Note that take complex conjugate since we require the inverse DFT
					DIM(IYY,IZZ)=DIM(IYY,IZZ)+
     +					COS(PHY+PHZ)*SQRT(AD(J))
					IM(IYY,IZZ)=IM(IYY,IZZ)-
     +					SIN(PHY+PHZ)*SQRT(AD(J))
				ENDIF
			    ENDIF
			ENDDO	
C Multiply array by checkered pattern of 1,-1 to shift origin in inverse
C DFT domain
			DO J=1,NDIFF
				IF(MOD(J,2).EQ.0) THEN
					JJ=-1
				ELSE
					JJ=1
				ENDIF
				DO K=1,NDIFF
					IF(MOD(K,2).EQ.0) THEN
						KK=-JJ
					ELSE
						KK=JJ
					ENDIF
					DIM(K,J)=DIM(K,J)*KK
					IM(K,J)=IM(K,J)*KK
				ENDDO
			ENDDO
C Take inverse DFT of principal plane distribution to find diffraction pattern
			NELS(1)=NDIFF
			NELS(2)=NDIFF
			NN=NELS(1)*NELS(2)
			IFAIL=-1
C			CALL C06FJF(2,NELS,NN,DIM,IM,WK,ND3,IFAIL)
			IF(IFAIL.NE.0) THEN
				ISTAT=-1
			ELSE
C Calculate power (intensity) of diffraction pattern
				DO J=1,NDIFF
					DO K=1,NDIFF
						DIM(K,J)=DIM(K,J)**2+IM(K,J)**2
					ENDDO
				ENDDO
			ENDIF
		ENDIF
	ELSE
		YBAR=0.0
		ZBAR=0.0
		RMS=0.0
	ENDIF
	END
