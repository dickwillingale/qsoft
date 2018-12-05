PRO QRI_ZOOMIMAGE, ima, xl, yl,zl

device,decomposed=0
; get QRI color table
;clut=qr_collut()

IF ISA(ima,'STRUCT') THEN xla = ima.xlab ELSE xla=""
IF ISA(ima,'STRUCT') THEN yla = ima.ylab ELSE yla="" 
IF ISA(ima,'STRUCT') THEN tit = ima.title ELSE tit="" 
IF ISA(ima,'STRUCT') THEN img = ima.data_array ELSE img = ima
IF ISA(ima,'STRUCT') THEN X_range = ima.xlim ELSE X_range = xl
IF ISA(ima,'STRUCT') THEN Y_range = ima.ylim ELSE Y_range = yl

  cgIMAGE,img, XRANGE = X_range, YRANGE = Y_range, XTITLE = xla, YTITLE = yla, $
    /axes, stretch = 1, minvalue = zl(0), maxvalue = zl(1), position = [0.15,0.15,0.85,0.85],$
    /normal,/reverse, CHARSIZE=2.5,title = tit, /KEEP_ASPECT_RATIO, color='black'

      cbarpos = [!X.Window[1], !Y.Window[0], !X.Window[1]+0.01, !Y.Window[1]]
  cgcolorbar,/vertical,/normal,range = zl,/right,position=cbarpos, CHARSIZE=1,color='black'

RETURN
END