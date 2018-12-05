# qR fitsio routines
if(!exists("qrf_init")) {
qrf_init<-function() {
	.Fortran("qrf_init")
	invisible()
}
qrf_init()
qr_fitsgetcol<-function(typ,ic,ii,ne) {
# Get column from table on fits file
	if(ne==0) return(NA)
	if(typ==1) {
		s<-.Fortran("qr_fitsgcvj",
		as.integer(NA_integer_),
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		vals=integer(length=ne),
		NAOK=T)
		s$vals
	} else if(typ==2) {
		s<-.Fortran("qr_fitsgcvd",
		as.double(NA_real_),
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		vals=double(length=ne),
		NAOK=T)
		s$vals
	} else if(typ==3) {	
		s<-.Fortran("qr_fitsgcvl",
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		vals=integer(length=ne))
		as.logical(s$vals)
	} else if(typ==4) {
		s<-.Fortran(
		"qr_fitsgcs",
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		val=character(length=1),
		is=integer(length=1),
		idf=integer(length=1))
		if(s$idf==1) {
			vals<-NA
		} else {
			vals<-substr(s$val,1,s$is)
		}
	} else if(typ==5) {
		nels2<-ne*2
		s<-.Fortran("qr_fitsgcvc",
		as.single(NA_real_),
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		vals=single(length=nels2),
		NAOK=T)
		dim(s$vals)<-c(2,ne)
		vals<-complex(real=s$vals[1,1:ne],imaginary=s$vals[2,1:ne])
	} else if(typ==6) {
		nels2<-ne*2
		s<-.Fortran("qr_fitsgcvm",
		as.double(NA_real_),
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		vals=double(length=nels2),
		NAOK=T)
		dim(s$vals)<-c(2,ne)
		vals<-complex(real=s$vals[1,1:ne],imaginary=s$vals[2,1:ne])
	} else if(typ==7) {
		s<-.Fortran("qr_fitsgcvj",
		as.integer(NA_integer_),
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		vals=integer(length=ne),
		NAOK=T)
		as.raw(s$vals)
	} else if(typ==8) {
		s<-.Fortran(
		"qr_fitsgcx",
		as.integer(ic),
		as.integer(ii),
		as.integer(ne),
		val=character(length=1))
		vals<-substr(s$val,1,ne)
	}
}
qr_fitshdu<- function(ihdu) {
# Set current fits HDU and return info
	a<-.Fortran("qr_fitshdu",
	as.integer(ihdu),
	as.integer(10),
	hdutype=integer(length=1),
	naxis=integer(length=1),
	naxes=integer(length=10),
	nkeys=integer(length=1))
	return(a)
}
qr_fitsopen<- function(filename) {
# Open fits file for read-only
	nn<-nchar(filename,type="chars")
	name<-substr(filename,1,nn)
	nrw<- 0
	a<-.Fortran("qr_fitsopen",
	as.integer(nn),
	as.character(name),
	as.integer(nrw),
	nhdu=integer(length=1))
	return(a)
}
qr_fitsread<-function(filename) {
# Read contents of fits file into R object
# Open file and get number of header units
	nhdu<-qr_fitsopen(filename)$nhdu
	dread<-date()
# Create container list
	if(nhdu==1) {
		fits<-list(FILENAME=filename,NHDU=nhdu,
		DATE=dread,FTRANS="qr_fitsread",
		primary=list())
	} else {
		fits<-list(FILENAME=filename,NHDU=nhdu,
		DATE=dread,FTRANS="qr_fitsread",
		primary=list(),extension=list())
	}
# Loop for all header units
	for(i in 1:nhdu) {
		j<-0
		nams<-c()
		coms<-vector(mode="character")
		icom<-0
		koms<-vector(mode="character")
		hdu<-list()
# Get info for this header unit
		a<-qr_fitshdu(i)
		hdutype<-a$hdutype
		naxis<-a$naxis
		nrows<-a$naxes[1]
		ncols<-a$naxes[2]
		if(hdutype==0) ncols<-1
		nkeys<-a$nkeys
# set hdutype
		j<-j+1
		hdu[[j]]<-hdutype
		nams[j]<-"HDUTYPE"
		koms[j]<- "0 image data, 1 ascii table, 2 binary table"
# Loop for all keywords
		for(k in 1:nkeys) {
			v<-.Fortran("qr_fitsgetkey",
			as.integer(k),
			key=character(length=1),
			ik=integer(length=1),
			sval=character(length=1),
			is=integer(length=1),
			jval=integer(length=1),
			dval=double(length=1),
			lval=logical(length=1),
			ktype=integer(length=1),
			cval=character(length=1),
			ic=integer(length=1))
			key<-substr(v$key,1,v$ik)
			if(v$ktype==1) {
				val<-v$jval
			} else if(v$ktype==2) {
				val<-v$dval
			} else if(v$ktype==3) {
				val<-v$lval
			} else if(v$ktype==4) {
				if(v$is>0) {
					val<-substr(v$sval,1,v$is)
				} else {
					val<-""
				}
			}
			if(key!=""&&key!="COMMENT"&&key!="HISTORY"
			&&key!="CONTINUE") {
				j<-j+1
				hdu[[j]]<-val
				nams[j]<-key
				if(v$ic>0) {
					koms[j]<-substr(v$cval,1,v$ic)
				} else {
					koms[j]<-""
				}
			} else {
				if(v$ic>0) {
			 		icom<-icom+1
					coms[icom]<-substr(v$cval,1,v$ic)
				}
				
			}
		}
# Get data types
		b<-.Fortran("qr_fitstypes",
		as.integer(hdutype),
		as.integer(ncols),
		ctype=integer(length=ncols),
		rp=integer(length=ncols))
		rp<-b$rp
		ctype<-b$ctype
		if(hdutype==0 && naxis>0) {
# Primary data array
			nels<-1
			for(ia in 1:naxis) nels<-nels*a$naxes[ia]
			if(ctype==1) {
				c<-.Fortran("qr_fitsgpvj",
				as.integer(NA_integer_),
				as.integer(nels),
				data_array=integer(length=nels),
				NAOK=T)
				j<-j+1
				hdu[[j]]<-array(c$data_array,
				dim=a$naxes[1:naxis])
				nams[j]<-"DATA_ARRAY"
			} else {
				c<-.Fortran("qr_fitsgpvd",
				as.double(NA_real_),
				as.integer(nels),
				data_array=double(length=nels),
				NAOK=T)
				j<-j+1
				hdu[[j]]<-array(c$data_array,
				dim=a$naxes[1:naxis])
				nams[j]<-"DATA_ARRAY"
			}
		} else if(naxis>0) {
# data table
			cnams<-list()
			ctab<-list()
			ict<-0
			rnams<-list()
			rtab<-list()
			irt<-0
			for(ic in 1:ncols) {
# get column name and variable repeat count if rp[ic]=0
				rr<-rp[ic]
				s<-.Fortran("qr_fitscolnam",
				as.integer(ic),
				as.integer(rr),
				as.integer(nrows),
				as.integer(255),
				colnam=character(length=1),
				iname=integer(length=1),
				vrep=integer(length=nrows))
				colnam<-substr(s$colnam,1,s$iname)
# check for FITS bit data
				if(ctype[ic]==8) {
					colnam=paste(colnam,"_bits",sep="")
				}
				vr<-s$vrep
# set parameters depending on column repeat and type
				if(rp[ic]==0) {
# Variable width column
					nc<-nrows
					ne<-vr
					cv<-list()
					qt<-4
				} else if(ctype[ic]==4 || ctype[ic]==8) {
# character or bit - get 1 R value per row
				        nc<-nrows
					ne<-rep(rp[ic],nrows)
					cv<-array(data=NA,dim=nrows)
					qt<-3
				} else if(rp[ic]>1) {
# Fixed width column >1
					nc<-1
					ne=rp[ic]*nrows
					cv=array(data=NA,
					dim=c(rp[ic],nrows))
					qt<-2
				} else {
# Column width 1
					nc<-1
					ne<-nrows
					qt<-1
				}
# Loop for calls to get columnn values
				for(ii in 1:nc) {
					vals<-qr_fitsgetcol(ctype[ic],
					ic,ii,ne[ii])
					if(length(vals)==0) {
						vals<-NA
					}
					if(qt==4) {
						cv[[ii]]=vals
					} else if(qt==3) {
						cv[ii]<-vals
					} else if(qt==2) {
						nn<-ne[ii]/nrows
						cv[1:nn,1:nrows]<-vals
					} else {
						cv<- vals
					}
				}
# Put values into correct table type
				if(qt==4 || qt==2) {
					irt<-irt+1
					rtab[[irt]]<-cv
					rnams[[irt]]<-colnam
				} else {
					ict<-ict+1
					ctab[[ict]]<-cv
					cnams[[ict]]<-colnam
				}
			}
# Put table into HDU
			if(ict>0) {
				j<-j+1
				names(ctab)<-cnams
				hdu[[j]]<-as.data.frame(ctab)
				nams[j]<-"table"
			}
			if(irt>0) {
				j<-j+1
				names(rtab)<-rnams
				hdu[[j]]<-rtab
				nams[j]<-"rtable"
			}
		}
# set list for comments
		j<-j+1
		nams[j]="comments"
		hdu[[j]]<-coms
		j<-j+1
		nams[j]="komments"
		hdu[[j]]<-koms
# Set names in HDU
		names(hdu)<-nams
# save in fits list
		if(i==1) {
			fits$primary=hdu
		} else {
			fits$extension[[i-1]]=hdu
		}
	}
	qr_fitsclose()
# return completed object
	return(fits)
}
qr_fitsprint<-function(fits) {
# Print/list contents of R object produced by qr_fitsread()
	maxels=100
	if(maxels>0) {
		maxsave<-getOption("max.print")
		options(max.print=maxels)
	}
	cat("filename       ",fits$FILENAME,"\n")
	cat("date           ",fits$DATE,"\n")
	cat("nhdu           ",fits$NHDU,"\n")
	cat("read by        ",fits$FTRANS,"\n")
	nhdu<-fits$NHDU
	for(i in 1:nhdu) {
		if(i==1) {
			hdu<-fits$primary
			cat("*** primary ***\n")
		} else {
			hdu<-fits$extension[[i-1]]
			cat("*** extension",i-1,"***\n")
		}
# List pure comment records
		coms<- hdu$comments
		ncom<- length(coms)
		if(ncom>0) {
			for(k in 1:ncom) cat(" ",coms[k],"\n")
		}
		koms<- hdu$komments
		nams<-names(hdu)
		for(j in 1:length(hdu)) {
			if(nams[j]=="DATA_ARRAY") {
				cat("DATA_ARRAY",typeof(hdu[[j]]),
				dim(hdu[[j]]),"\n")
				DATA_ARRAY<-hdu$DATA_ARRAY
				cat(DATA_ARRAY[1:1,1:1],"...\n")
				cat("minval",min(DATA_ARRAY),
				"maxval",max(DATA_ARRAY),"\n")
			} else if(nams[j]=="table") {
				at<-hdu[[j]]
				nc<-length(at)
				cat("table",nc,"columns\n")
				print(at)
			} else if(nams[j]=="rtable") {
				rt<-hdu[[j]]
				inew<-0
				nc<-length(rt)
				cat("rtable",nc,"columns\n")
				cnams<-names(rt)
				for(k in 1:nc) {
					newrt<-list()
					newnams<-list()
					col<-rt[[k]]
					nd<-dim(col)
					rp<-nd[1]
					nrows<-nd[2]
					cat(cnams[k],rp,nrows,"\n")
					print(col[1:min(rp,10),1:min(nrows,10)])
				}
			} else if(nams[j]!="komments"&&nams[j]!="comments") {
				if(koms[j]!="") {
					cat(nams[j],"=",hdu[[j]]," : ",
					koms[j],"\n")
				} else {
					cat(nams[j],"=",hdu[[j]],"\n")
				}
			}
		}
	}
	if(maxels>0) options(max.print=maxsave)
}
qr_fitsnew<-function(fname) {
# Create new fits file
	.Fortran("qr_fitsnew",
	as.integer(nchar(fname)),
	as.character(fname))
	invisible()
}
qr_fitssave<-function(fits,fname) {
# Save as fits file the contents of an R object with same format as
# produced by qr_fitsread()
	qr_fitsnew(fname)
	fits$FILENAME<- fname
	fits$DATE<- date()
	iext<- grep("extension",names(fits))
	if(length(iext)==0) iext<- 0
	if(iext>0) {
		nhdu<- length(fits$extension)+1
	} else {
		nhdu<- 1
	}
	fits$NHDU<- nhdu
	fits$FTRANS<- "qr_fitssave()"
# Loop for all HDU
	for(i in 1:nhdu) {
		if(i==1) {
			hdu<-fits$primary
		} else {
			hdu<-fits$extension[[i-1]]
		}
		nams<- names(hdu)
		idata<-grep("DATA_ARRAY",nams)
		if(idata==0) idata<- 0
		if(idata>0) {
			obj<- hdu[[idata]]
			typ <- typeof(obj)
			if(typ=="logical") {
                        	la<- as.integer(obj)
                        	dim(la)<- dim(obj)
                        	qr_fitsparrj(la)
               		 } else if(typ=="integer") {
                        	qr_fitsparrj(obj)
                	} else if(typ=="double") {
                       		qr_fitsparrd(obj)
               		}
		}
		itab<- grep("table",nams)
		if(length(itab)==0) itab<- 0
		irtab<<- grep("rtable",nams)
		if(length(irtab)==0) irtab<- 0
		if(itab==0) {
			itab<- irtab
			irtab<- 0
 		}
		if(itab>0) {
			obj<- hdu[[itab]]
			len<- length(obj)
			nd<- dim(obj)
			nam<- names(obj)
# find data type of each column in table
			ity<- numeric(length=nd[2])
			ire<- numeric(length=nd[2])
			iwi<- numeric(length=nd[2])
			for(i in 1:len) {
			    if(typeof(obj[[i]])=="integer") {
				ity[i]<- 1
				ire[i]<- 1
				iwi[i]<- 4
			    } else if(typeof(obj[[i]])=="double") {
				ity[i]<- 2
				ire[i]<- 1
				iwi[i]<- 8
			    } else if(typeof(obj[[i]])=="logical") {
				ity[i]<- 3
				ire[i]<- 1
				iwi[i]<- 4
			    } else if(typeof(obj[[i]])=="character") {
				if(grepl("_bits",nam[i])) {
					ity[i]<- 8
					ire[i]<- nchar(max(obj[[i]]))
					iwi[i]<- 0
				} else {
					ity[i]<- 4
					ire[i]<- 1
					iwi[i]<- 0
				}
				iwi[i]<- max(nchar(obj[[i]]))
			    } else if(typeof(obj[[i]])=="complex") {
				ity[i]<- 5
				ire[i]<- 1
				iwi[i]<- 8
			    } else if(typeof(obj[[i]])=="raw") {
				ity[i]<- 7
				ire[i]<- 1
				iwi[i]<- 1
			    }
			}
		}
		if(irtab>0) {
			objr<- hdu[[irtab]]
			lenr<- length(objr)
			namr<- names(objr)
			nd[2]<- nd[2]+lenr
# find data type of each column in rtable
			for(ii in 1:lenr) {
			    nd[1]<- max(nd[1],length(objr[[ii]]))
			    i<- len+ii
			    if(typeof(objr[[ii]])=="integer") {
				ity[i]<- 1
				ire[i]<- 1
				iwi[i]<- 4
			    } else if(typeof(objr[[ii]])=="double") {
				ity[i]<- 2
				ire[i]<- 1
				iwi[i]<- 8
			    } else if(typeof(objr[[ii]])=="logical") {
				ity[i]<- 3
				ire[i]<- 1
				iwi[i]<- 4
			    } else if(typeof(objr[[ii]])=="character") {
				if(grepl("_bits",namr[ii])) {
					ity[i]<- 8
					ire[i]<- nchar(max(objr[[ii]]))
					iwi[i]<- 0
				} else {
					ity[i]<- 4
					ire[i]<- 1
					iwi[i]<- 0
				}
				iwi[i]<- max(nchar(objr[[ii]]))
			    } else if(typeof(objr[[ii]])=="complex") {
				ity[i]<- 5
				ire[i]<- 1
				iwi[i]<- 8
			    } else if(typeof(objr[[ii]])=="raw") {
				ity[i]<- 7
				ire[i]<- 1
				iwi[i]<- 1
			    }
			}
# Set binary table header in FITS file
			qr_fitspbtab(nd[1],nd[2],ity,ire,iwi)
		}
# Put data into binary table a column at a time
		if(itab>0) {
		    for(i in 1:len) {
			cat("column",i,nam[[i]],ity[i],"\n")
			if(ity[i]==1) {
				qr_fitspcolj(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==2) {
				qr_fitspcold(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==3) {
				qr_fitspcoll(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==4) {
				nn<- length(obj[[i]])
				for(ir in 1:nn) {
				  qr_fitspcols(i,ir,obj[[i]][ir],nam[[i]],"")
				}
			} else if(ity[i]==5) {
				qr_fitspcolc(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==7) {
				qr_fitspcolb(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==8) {
				nn<- length(obj[[i]])
				fnm<- substr(nam[[i]],1,nchar(nam[[i]])-5)
				for(ir in 1:nn) {
				  qr_fitspcolx(i,ir,obj[[i]][ir],fnm,"")
				}
                        }
		    }
		}
		if(irtab>0) {
		    for(ii in 1:lenr) {
			i<- len+ii
			if(ity[i]==1) {
				qr_fitspcolj(i,objr[[ii]],namr[[ii]],"")
			} else if(ity[i]==2) {
				qr_fitspcold(i,objr[[ii]],namr[[ii]],"")
			} else if(ity[i]==3) {
				qr_fitspcoll(i,objr[[ii]],namr[[ii]],"")
			} else if(ity[i]==4) {
				nn<- length(objr[[ii]])
				for(ir in 1:nn) {
				  qr_fitspcols(i,ir,objr[[ii]][ir],namr[[ii]],"")
				}
			} else if(ity[i]==5) {
				qr_fitspcolc(i,objr[[ii]],namr[[i]],"")
			} else if(ity[i]==7) {
				qr_fitspcolb(i,objr[[ii]],namr[[ii]],"")
			} else if(ity[i]==8) {
				nn<- length(objr[[ii]])
				fnm<- substr(namr[[ii]],1,nchar(namr[[ii]])-5)
				for(ir in 1:nn) {
				  qr_fitspcolx(i,ir,objr[[ii]][ir],fnm,"")
				}
                        }
		    }
		}
# Save pure comment records
		coms<- hdu$comments
		ncom<- length(coms)
		if(ncom>0) {
			for(k in 1:ncom) {
				qr_fitspcom(coms[k])
			}
		}
# Loop to pick out keywords
		for(j in 1:length(hdu)) {
		    if(nams[j]!="DATA_ARRAY" && nams[j]!="table"
		    && nams[j]!="rtable" && nams[j]!="comments"
		    && nams[j]!="kcomments") {
			obj<- hdu[[j]]
			typ<- typeof(obj)
			objnam<- nams[j]
			com<- hdu$Kcomments[j]
			if(typ=="logical") {
				qr_fitspkeyl(objnam,obj,com)
			} else if(typ=="integer") {
				qr_fitspkeyj(objnam,obj,com)
			} else if(typ=="double") {
				qr_fitspkeyd(objnam,obj,com)
			} else if(typ=="complex") {
				rpart<- paste(objnam,"_realpart",sep="")
				qr_fitspkeyd(rpart,Re(obj),com)
				ipart<- paste(objnam,"_imagpart",sep="")
				qr_fitspkeyd(ipart,Im(obj),com)
			} else if(typ=="raw") {
				qr_fitspkeyj(objnam,as.integer(obj),com)
			} else if(typ=="character") {
				qr_fitspkeys(objnam,obj,com)
			} else {
				qr_fitspkeys(objnam,typ,com)
			}
		    }
		}
	}
	qr_fitsclose()
}
qr_fitsempty<-function() {
# Write empty HDU to fits file
	.Fortran("qr_fitsempty")
	invisible()
}
qr_fitsparrd<-function(arr) {
# Put double array onto fits file as primary or extension data array
	nels<-dim(arr)
	ndim<-length(nels)
	nel<-length(arr)
	.Fortran("qr_fitsparrd",
	as.integer(ndim),
	as.integer(nels),
	as.integer(nel),
	as.double(arr))
	invisible()
}
qr_fitsparrj<-function(iarr) {
# Put integer array onto fits file as primary or extension data array
	nels<-dim(iarr)
	ndim<-length(nels)
	nel<-length(iarr)
	.Fortran("qr_fitsparrj",
	as.integer(ndim),
	as.integer(nels),
	as.integer(nel),
	as.integer(iarr))
	invisible()
}
qr_fitspbtab<-function(nrows,ncols,ity,ire,iwi) {
# Create binary table HDU on fits file
	.Fortran("qr_fitspbtab",
	as.integer(nrows),
	as.integer(ncols),
	as.integer(ity),
	as.integer(ire),
	as.integer(iwi))
	invisible()
}
qr_fitspcolb<-function(ncol,col,cname,cunit) {
# Put raw byte column data on to fits file
	.Fortran("qr_fitspcolb",
	as.integer(ncol),
	as.integer(length(col)),
	as.integer(col),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}	
qr_fitspcolc<-function(ncol,col,cname,cunit) {
# Put complex column data on to fits file
	len<- length(col)
	carr<- array(dim=c(2,len)) 
	carr[1,1:len]<- Re(col)
	carr[2,1:len]<- Im(col)
	.Fortran("qr_fitspcolc",
	as.integer(ncol),
	as.integer(len),
	as.single(carr),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}	
qr_fitspcold<-function(ncol,col,cname,cunit) {
# Put double precision column data on to fits file
	.Fortran("qr_fitspcold",
	as.integer(ncol),
	as.integer(length(col)),
	as.double(col),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}	
qr_fitspcolj<-function(ncol,col,cname,cunit) {
# Put integer column data on to fits file
	.Fortran("qr_fitspcolj",
	as.integer(ncol),
	as.integer(length(col)),
	as.integer(col),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}
qr_fitspcoll<-function(ncol,col,cname,cunit) {
# Put logical column data on to fits file
	.Fortran("qr_fitspcoll",
	as.integer(ncol),
	as.integer(length(col)),
	as.integer(col),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}
qr_fitspcols<-function(ncol,nrow,col,cname,cunit) {
# Put single row string of column data on to fits file
	.Fortran("qr_fitspcols",
	as.integer(ncol),
	as.integer(nrow),
	as.integer(nchar(col)),
	as.character(col),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}
qr_fitspcolx<-function(ncol,nrow,col,cname,cunit) {
# Put single row of bits to fits file
	.Fortran("qr_fitspcolx",
	as.integer(ncol),
	as.integer(nrow),
	as.integer(nchar(col)),
	as.character(col),
	as.integer(nchar(cname)),
	as.character(cname),
	as.integer(nchar(cunit)),
	as.character(cunit))
	invisible()
}
qr_fitspcom<-function(com) {
# Put comment card on to fits file
	.Fortran("qr_fitspcom",
	as.integer(nchar(com)),
	as.character(com))
	invisible()
}
qr_fitsphis<-function(his) {
# Put history card on to fits file
	.Fortran("qr_fitsphis",
	as.integer(nchar(his)),
	as.character(his))
	invisible()
}
qr_fitspkeys<-function(key,str,com) {
# Put string keyword on to fits file
	.Fortran("qr_fitspkeys",
	as.character(key),
	as.character(str),
	as.character(com))
	invisible()
}
qr_fitspkeyd<-function(key,val,com) {
# Put double precision keyword on to fits file
	.Fortran("qr_fitspkeyd",
	as.character(key),
	as.double(val),
	as.character(com))
	invisible()
}
qr_fitspkeyj<-function(key,ival,com) {
# Put integer keyword on to fits file
	.Fortran("qr_fitspkeyj",
	as.character(key),
	as.integer(ival),
	as.character(com))
	invisible()
}
qr_fitspkeyl<-function(key,lval,com) {
# Put logical keyword on to fits file
	.Fortran("qr_fitspkeyl",
	as.character(key),
	as.logical(lval),
	as.character(com))
	invisible()
}
qr_fitspobj<- function(obj,objnam) {
# Put elements of R object structure into fits file
	typ<- typeof(obj)
	cla<- class(obj)
	len<- length(obj)
	nam<- names(obj)
	if(cla=="data.frame") {
# R data frame entered as a FITS binary table extension
		nd<- dim(obj)
		cat("R data frame",objnam,"nrows:",nd[1],"ncols:",nd[2],"\n")
# find data type of each column in table
		ity<- numeric(length=nd[2])
		ire<- numeric(length=nd[2])
		iwi<- numeric(length=nd[2])
		for(i in 1:len) {
			if(typeof(obj[[i]])=="integer") {
				ity[i]<- 1
				ire[i]<- 1
				iwi[i]<- 4
			} else if(typeof(obj[[i]])=="double") {
				ity[i]<- 2
				ire[i]<- 1
				iwi[i]<- 8
			} else if(typeof(obj[[i]])=="logical") {
				ity[i]<- 3
				ire[i]<- 1
				iwi[i]<- 4
			} else if(typeof(obj[[i]])=="character") {
				if(grepl("_bits",nam[i])) {
					ity[i]<- 8
					ire[i]<- nchar(max(obj[[i]]))
					iwi[i]<- 0
				} else {
					ity[i]<- 4
					ire[i]<- 1
					iwi[i]<- 0
				}
				iwi[i]<- max(nchar(obj[[i]]))
			} else if(typeof(obj[[i]])=="complex") {
				ity[i]<- 5
				ire[i]<- 1
				iwi[i]<- 8
			} else if(typeof(obj[[i]])=="raw") {
				ity[i]<- 7
				ire[i]<- 1
				iwi[i]<- 1
			}
		}
# Check for row.names
		rnames<- row.names(obj)
		lrn<- length(rnames)
		if(lrn==nd[1]) {
			nd[2]<- nd[2]+1
			ity[len+1]<- 4
			ire[len+1]<- 1
			iwi[len+1]<- max(nchar(rnames))
		}
# Set binary table header in FITS file
		qr_fitspbtab(nd[1],nd[2],ity,ire,iwi)
		com<- paste("R class",cla)
		if(length(objnam)>0) qr_fitspkeys("EXTNAME",objnam,com)
# Put data into binary table a column at a time
		nam<- names(obj)
		for(i in 1:len) {
			cat("column",i,nam[[i]],ity[i],"\n")
			if(ity[i]==1) {
				qr_fitspcolj(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==2) {
				qr_fitspcold(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==3) {
				qr_fitspcoll(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==4) {
				nn<- length(obj[[i]])
				for(ir in 1:nn) {
				  qr_fitspcols(i,ir,obj[[i]][ir],nam[[i]],"")
				}
			} else if(ity[i]==5) {
				qr_fitspcolc(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==7) {
				qr_fitspcolb(i,obj[[i]],nam[[i]],"")
			} else if(ity[i]==8) {
				nn<- length(obj[[i]])
				fnm<- substr(nam[[i]],1,nchar(nam[[i]])-5)
				for(ir in 1:nn) {
				  qr_fitspcolx(i,ir,obj[[i]][ir],fnm,"")
				}
                        }
		}
# Put row.names column if required
		if(lrn==nd[1]) {
			for(ir in 1:lrn) {
			  qr_fitspcols(nd[2],ir,rnames[ir],"row.names","")
			}
		}	
	} else if(cla=="list") {
# R list - create an empty extension and save name as EXTNAME
		cat("R list",objnam,"\n")
		qr_fitsempty()
		com<- paste("R class",cla)
		if(length(objnam)>0) qr_fitspkeys("EXTNAME",objnam,com)
# Work through list
		nam<- names(obj)
		for(i in 1:len) {
                        qr_fitspobj(obj[[i]],nam[[i]])
                }
	} else if(len==1) {
# scalar primitive entered as a FITS keyword
		cat("keyword",objnam,typ,"\n")
		com<- paste("R type",typ)
		if(typ=="logical") {
			qr_fitspkeyl(objnam,obj,com)
		} else if(typ=="integer") {
			qr_fitspkeyj(objnam,obj,com)
		} else if(typ=="double") {
			qr_fitspkeyd(objnam,obj,com)
		} else if(typ=="complex") {
			rpart<- paste(objnam,"_realpart",sep="")
			qr_fitspkeyd(rpart,Re(obj),com)
			ipart<- paste(objnam,"_imagpart",sep="")
			qr_fitspkeyd(ipart,Im(obj),com)
		} else if(typ=="raw") {
			qr_fitspkeyj(objnam,as.integer(obj),com)
		} else if(typ=="character") {
			qr_fitspkeys(objnam,obj,com)
		} else {
			qr_fitspkeys(objnam,typ,com)
		}
	} else {
# primitive array entered as primary FITS array or FITS extension array
		cat("data array",objnam,typ,"\n")
		if(typ=="logical") {
			la<- as.integer(obj)
			dim(la)<- dim(obj)
			qr_fitsparrj(la)
		} else if(typ=="integer") {
			qr_fitsparrj(obj)
		} else if(typ=="double") {
			qr_fitsparrd(obj)
		}
		com<- paste("R class",cla)
		if(length(objnam)>0) qr_fitspkeys("EXTNAME",objnam,com)
	}
	invisible()
}
qr_fitswrite<-function(obj,fname) {
# Write R object structure to new fits file
	qr_fitsnew(fname)
	objnam<- deparse(substitute(obj))
	qr_fitspobj(obj,objnam)
	qr_fitsclose()
	invisible()
}
qr_fitsupdate<-function(filename) {
# Open fits file and get number of header units
	nn<-nchar(filename,type="chars")
	name<-substr(filename,1,nn)
	nrw<- 1
	nhdu<-.Fortran("qr_fitsopen",as.integer(nn),as.character(name),
	as.integer(nrw),nhdu=integer(length=1))$nhdu
# Return number of Header Data Units
	return(nhdu)
}
qr_fitsclose<- function() {
	.Fortran("qr_fitsclose")
	invisible()
}
# Quick get out of jail
  qq<-function() {
    q(save="no")
    invisible()
  }
}
