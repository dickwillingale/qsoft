*+SRT_LISTMPOARR    List MPO array parameters 
        SUBROUTINE SRT_LISTMPOARR(IU)
        IMPLICIT NONE
        INTEGER IU
*IU  input    output channel
*-Author Dick Willingale 2019-Mar-6
        INCLUDE 'SRT_COM'
        INTEGER I,II,ILINE,J
        PARAMETER (ILINE=4)
999     FORMAT(28X,52('-'))
1000    FORMAT(9X,I4,' MPOs in array')
1001    FORMAT(22X,4(1X,I12))
1002    FORMAT(9X,'x position         ',4(1X,G12.5))
1003    FORMAT(9X,'y position         ',4(1X,G12.5))
1004    FORMAT(9X,'width              ',4(1X,G12.5))
1005    FORMAT(9X,'height             ',4(1X,G12.5))
1006    FORMAT(9X,'rotation angle     ',4(1X,G12.5))
1007    FORMAT(9X,'thickness          ',4(1X,G12.5))
1008    FORMAT(9X,'radius of curvature',4(1X,G12.5))
1009    FORMAT(9X,'multifibre size    ',4(1X,G12.5))
1010    FORMAT(9X,'channel pitch      ',4(1X,G12.5))
1011    FORMAT(9X,'wall thickness     ',4(1X,G12.5))
1012    FORMAT(9X,'surface quality    ',4(1X,G12.5))
1013    FORMAT(9X,'bias angle error x ',4(1X,G12.5))
1014    FORMAT(9X,'bias angle error y ',4(1X,G12.5))
1015    FORMAT(9X,'efficiency         ',4(1X,G12.5))
1016    FORMAT(9X,'spare              ',4(1X,G12.5))
C Number of MPOs
        WRITE(IU,*)
        WRITE(IU,1000) nrofelts
# Loop for all MPOs
        DO I=1,nrofelts,ILINE
                II=MIN(I+ILINE-1,nrofelts)
                WRITE(IU,999)
                WRITE(IU,1001) (J,J=I,II)
                WRITE(IU,999)
C Position of module centre
                WRITE(IU,1002) (rlist(J),J=I,II)
                WRITE(IU,1003) (philist(J),J=I,II)
C Dimensions of module aperture
                WRITE(IU,1004) (wlist(J),J=I,II)
                WRITE(IU,1005) (hlist(J),J=I,II)
C Rotation angle
                WRITE(IU,1006) (thlist(J),J=I,II)
C Thickness
                WRITE(IU,1007) (llist(J),J=I,II)
C RoC
                WRITE(IU,1008) (clist(J),J=I,II)
C Multifibre size
                WRITE(IU,1009) (glist(J),J=I,II)
C Channel pitch
                WRITE(IU,1010) (olist(J),J=I,II)
C Wall thickness
                WRITE(IU,1011) (plist(J),J=I,II)
C Surface quality
                WRITE(IU,1012) (qlist(J),J=I,II)
C Bias angles
                WRITE(IU,1013) (ulist(J),J=I,II)
                WRITE(IU,1014) (vlist(J),J=I,II)
C Efficiency wrt theory
                WRITE(IU,1015) (zlist(J),J=I,II)
C Spare parameter
                WRITE(IU,1016) (slist(J),J=I,II)
        ENDDO
        WRITE(IU,*)
        END
