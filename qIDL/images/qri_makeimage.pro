FUNCTION QRI_MAKEIMAGE, data_array,xleft,xright,ybot,ytop

  nx = LONG(N_ELEMENTS(data_array(*,0)))
  ny = LONG(N_ELEMENTS(data_array(0,*)))


 B = QRI_SETFIELD(nx,DOUBLE(xleft),DOUBLE(xright),ny,DOUBLE(ybot),DOUBLE(ytop))

  
C = {data_array:data_array, $
  nx:nx,xleft:xleft,xright:xright,xsam:b.xsam,$
  ny:ny,ybot:ybot,ytop:ytop,ysam:b.ysam,$
  xp:b.xp,yp:b.yp,xlim:[xleft,xright],ylim:[ybot,ytop],$
  zlim:[MIN(data_array),MAX(data_array)],xlab:'',ylab:'',title:''}

RETURN, C
END