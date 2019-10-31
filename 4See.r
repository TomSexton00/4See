#############################################################
#4See: A GUI for dynamic browsing of 4C data
#Copyright (C) (2019) Tom Sexton & Yousra Ben Zouari
#
#This program is free software; you can distribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
##############################################################


options(warn=-1)
suppressPackageStartupMessages(require(tcltk2))
suppressPackageStartupMessages(require(tkrplot))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(caTools))
suppressPackageStartupMessages(require(rtracklayer))

#Note of called 'global' variables
#Viewpoint information
cis.chrom=NULL
vp.pos=NULL
#settings - list for each condition level: name (ID dataset name), entries (all names for level), color
settings=NULL
#data - the list of 4C cis files
data=NULL
#set - user input for managing replicates, used to generate settings
set=list(lab=list(),entry=list())
#pset - user input for managing conditions for plotting
pset=list()
#ints - interaction sets
ints=list()
#pints - tmp store of interaction parameters
pint=list(name=list(),plot=list())
#intplots - tmp store of if interaction is plotted or not
intplots=list()
#genes track
genome=NULL
#epigenomic tracks
epigenome=list()
#trks - tmp store of track parameters
trks=list(lab=list(),entry=list())
#cols, the default colors for each plot
plotcols=c("black","red","blue","green","cyan","purple","yellow","gray")


###################################################################################
#Called interaction overlays

#Change color helper - interactions
changeintcolor=function() {
	to.change = as.integer(tclvalue(tkget(to.change)))
	color = tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=ints[[to.change]]$color, title="Choose Color"))))
	tkconfigure(manint$env$canvas[[to.change]], bg=color)
	if(nchar(color)>0) {
		ints[[to.change]]$color <<- color
		tkconfigure(manint$env$canvas[[to.change]], bg=color)
	}
}

#Interaction plot settings util
verify.int.settings=function() {
	for (i in 1:length(ints)) {
		ints[[i]]$name <<- as.character(tkget(pint$name[[i]]))
		ints[[i]]$plot <<- as.character(tclvalue(intplots[[i]]))
	}
	tkdestroy(manint)
}

#Interaction plot settings
manage.ints=function() {
	manint <<- tktoplevel()
	tktitle(manint)="Plot settings for interactions"
	labcol=tklabel(manint,text="Color")
	labplot=tklabel(manint,text="Plot")
	manint$env$change = tk2button(manint, text="Change Color", command=changeintcolor)
	to.change <<- tkentry(manint,width=5,textvariable=tclVar("1"))
	tkgrid(labcol, row=1, column=2, padx=15, pady=5)
	tkgrid(labplot, row=1, column=3, padx=15, pady=5)
	tkgrid(manint$env$change, row=1, column=4, padx=15, pady=5)
	tkgrid(to.change, row=1, column=5, padx=15, pady=5)
	for (i in 1:length(ints)) {
		pint$name[[i]] <<- tkentry(manint,width=10,textvariable=tclVar(ints[[i]]$name))
		manint$env$canvas[[i]] = tk2canvas(manint,width=30,height=10,bg=ints[[i]]$color)
		manint$env$plot[[i]] <<- tk2checkbutton(manint)
		intplots[[i]] <<- tclVar(ints[[i]]$plot)
		tkconfigure(manint$env$plot[[i]], variable=intplots[[i]])
		tkgrid(pint$name[[i]],row=i+1,column=1,padx=15,pady=5)
		tkgrid(manint$env$canvas[[i]],row=i+1,column=2,padx=15,pady=5)
		tkgrid(manint$env$plot[[i]],row=i+1,column=3,pady=5)
	}
	done=tkbutton(manint,text="OK",command=verify.int.settings)
	tkgrid(done,row=length(ints)+2,column=1,pady=5)
}

#Load interactions
load.ints=function() {
	files=tclvalue(tkgetOpenFile(filetypes = "{{Peak files} {.bed}} {{All files} *}",multiple=TRUE))
	if(!nchar(files)) {
		return()
	}
	files=unlist(strsplit(files,split=" ",fixed=FALSE,perl=FALSE,useBytes=FALSE))
	for (i in 1:length(files)) {
		ints[[i]] <<- list()
		int.name = unlist(strsplit(files[i],"/"))
		int.name=int.name[length(int.name)]
		int.name=unlist(strsplit(int.name,"\\."))
		int.name=paste(int.name[-(length(int.name))],collapse="")
		ints[[i]]$name <<- int.name
		ints[[i]]$color <<- plotcols[i]
		ints[[i]]$plot <<- "1"
		ints[[i]]$bed <<- read.table(files[i],header=TRUE,stringsAsFactors=FALSE)
	}
}

##################################################################
#Genes and epigenomic tracks

#Gene track loading
load.genes=function() {
	fn = tclvalue(tkgetOpenFile(filetypes = "{{gene files} {.txt}} {{all files} *}"))
	if(!nchar(fn)) {
		return()
	}
	genome <<- read.table(fn,header=TRUE,stringsAsFactors=FALSE)
}

#Epigenomic tracks

#Change color helper function - tracks
changetrkcolor=function() {
	to.change = as.integer(tclvalue(tkget(to.change)))
	color = tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=epigenome[[to.change]]$color, title="Choose Color"))))
	tkconfigure(mantrk$env$canvas[[to.change]], bg=color)
	if(nchar(color)>0) {
		epigenome[[to.change]]$color <<- color
		tkconfigure(mantrk$env$canvas[[to.change]], bg=color)
	}
}

#Track plot setting util
verify.trk.settings=function() {
	for (i in 1:length(epigenome)) {
		level = as.integer(tkget(trks$entry[[i]]))
		epigenome[[i]]$level <<- level
	}
	tkdestroy(mantrk)
}

#Plotting setting for tracks
manage.tracks=function() {
	mantrk <<- tktoplevel()
	tktitle(mantrk)="Plot settings for tracks"
	labcol=tklabel(mantrk,text="Color")
	lablev=tklabel(mantrk,text="Level")
	mantrk$env$change = tk2button(mantrk,text="Change Color", command=changetrkcolor)
	to.change <<- tkentry(mantrk,width=5,textvariable=tclVar("1"))
	tkgrid(labcol,row=1,column=2,padx=15,pady=5)
	tkgrid(lablev,row=1,column=3,padx=15,pady=5)
	tkgrid(mantrk$env$change,row=1,column=4,padx=15,pady=5)
	tkgrid(to.change,row=1,column=5,padx=15,pady=5)
	for (i in 1:length(epigenome)) {
		trks$lab[[i]] <<- tklabel(mantrk,text=names(epigenome)[i])
		mantrk$env$canvas[[i]] = tk2canvas(mantrk,width=30,height=10,bg=epigenome[[i]]$color)
		trks$entry[[i]] <<- tkentry(mantrk,width=5,textvariable=tclVar(epigenome[[i]]$level))
		tkgrid(trks$lab[[i]],row=i+1,column=1,padx=15,pady=5)
		tkgrid(mantrk$env$canvas[[i]],row=i+1,column=2,padx=15,pady=5)
		tkgrid(trks$entry[[i]],row=i+1,column=3,padx=15,pady=5)
	}
	done=tkbutton(mantrk,text="OK",command=verify.trk.settings)
	tkgrid(done,row=length(epigenome)+2,column=1,pady=5)
}

#Import bigwig/bedgraph tracks
import.bw=function(files) {
	for (i in 1:length(files)) {
		cat(paste("Reading file:",files[i],"\n"))
		trackname = unlist(strsplit(files[i],split="/",fixed=FALSE,perl=FALSE,useBytes=FALSE))
		trackname = trackname[length(trackname)]
		trackname = unlist(strsplit(trackname,split="\\."))[1]
		epigenome[[trackname]] <<- list()
		epigenome[[trackname]]$level <<- 0
		epigenome[[trackname]]$color <<- plotcols[i]
		epigenome[[trackname]]$track <<- import(files[i])
		cat("Finished reading file\n")
	}
}

#Epigenomic track load wrapper
load.epigenome=function() {
	input = tclvalue(tkgetOpenFile(filetypes = "{{Track Files} {.bw .bedGraph}} {{All files} *}",multiple=TRUE))
	if(!nchar(input)) {
		return()
	}
	files = unlist(strsplit(input,split=" ",fixed=FALSE,perl=FALSE,useBytes=FALSE))
	import.bw(files)
	manage.tracks()
}




###################################################################################
#4C profile inputs


#Change color helper function - condition/replicate set
changecolor=function() {
	to.change = as.integer(tclvalue(tkget(to.change)))
	color = tclvalue(.Tcl(paste("tk_chooseColor",.Tcl.args(initialcolor=settings[[to.change]]$color, title="Choose Color"))))
	tkconfigure(manage$env$canvas[[to.change]], bg=color)
	if(nchar(color)>0) {
		settings[[to.change]]$color <<- color
		tkconfigure(manage$env$canvas[[to.change]], bg=color)
	}
}


#Condition/replicate plot setting util
verify.plot.settings=function() {
	for (i in 1:length(settings)) {
		settings[[i]]$name <<- as.character(tkget(pset[[i]]))
	}
	tkdestroy(manage)
}

#Plotting settings for condition/replicates
manage.set=function() {
	manage <<- tktoplevel()
	tktitle(manage)="Plot settings for conditions"
	labname=tklabel(manage,text="Name")
	labcol=tklabel(manage,text="Color")
	manage$env$change = tk2button(manage, text="Change Color", command = changecolor)
	to.change <<- tkentry(manage,width=5,textvariable=tclVar("1"))
	tkgrid(labname, row=1, column=2, padx=15, pady=5)
	tkgrid(labcol, row=1, column=3, padx=15, pady=5)
	tkgrid(manage$env$change, row=1, column=4, padx = 15, pady=5)
	tkgrid(to.change, row=1, column=5, padx=5, pady=5)
	for (i in 1:length(settings)) {
		lablevel=tklabel(manage, text=paste(settings[[i]][["entries"]],collapse=","))
		tkgrid(lablevel, row=i+1, column=1, sticky="e",pady=5)
		pset[[i]] <<- tkentry(manage,width=10,textvariable=tclVar(settings[[i]]$name))
		manage$env$canvas[[i]] = tk2canvas(manage, width=30, height=10, bg=settings[[i]]$color)
		tkgrid(pset[[i]],row=i+1, column=2, padx=15, pady=5)
		tkgrid(manage$env$canvas[[i]],row=i+1,column=3, padx=15, pady=5)
	}
	done=tkbutton(manage,text="OK",command=verify.plot.settings)
	tkgrid(done,row=length(settings)+2,column=1,pady=5)
}


#Condition/replicate setting util
verify.settings=function() {
	assign=NULL
	indices=NULL
	for (i in 1:length(data)) {
		level = as.integer(tkget(set$entry[[i]]))
		assign[i] = level
	}
	names(assign) = names(data)
	assignb=assign[is.integer(assign) & assign>0]
	levels=unique(assignb)
	settings <<- list()
	for (l in 1:length(levels)) {
		settings[[levels[l]]] <<- list()
		hits = which(assign==levels[l])
		settings[[levels[l]]][["name"]] <<- names(hits)[1]
		settings[[levels[l]]][["indices"]] <<- hits
		settings[[levels[l]]][["entries"]] <<- names(hits)
		settings[[levels[l]]][["color"]] <<- plotcols[l]
	}
	tkdestroy(select)
	return(settings)
}


#Allocate 4C sets to condition/replicates
set.conditions=function(data) {
	select <<- tktoplevel()
	tktitle(select)="Manage conditions and replicates"
	for (i in 1:length(data)) {
		set$lab[[i]] <<- tklabel(select,text=names(data)[i])
		set$entry[[i]] <<- tkentry(select,width=10,textvariable=tclVar("1"))
		tkgrid(set$lab[[i]],row=i,column=1,sticky="e",pady=5)
		tkgrid(set$entry[[i]],row=i,column=2,sticky="w",pady=5)
	}
	done=tkbutton(select,text="OK",command=verify.settings)
	tkgrid(done,row=length(data)+1,column=1,pady=5)
}

#Open 4C cis file(s) from dialog box, if multiple, check for same vp and auto prompt for conditions setting
load.4c = function() {
	ifiles = tclvalue(tkgetOpenFile(filetypes = "{{4C files} {.cis}} {{all files} *}",multiple=TRUE))
	if(!nchar(ifiles)) {
		return()
	}
	data <<- NULL
	ifiles=unlist(strsplit(ifiles,split=" ",fixed=FALSE,perl=FALSE,useBytes=FALSE))
	if (length(ifiles)==1) {
		con=file(ifiles[1],"r")
		info=readLines(con,n=1)
		close(con)
		info=unlist(strsplit(info,"\\t"))
		cis.chrom <<- info[2]
		vp.pos <<- as.numeric(info[3])
		tab=read.table(ifiles[1],header=FALSE,skip=1)
		colnames(tab)=c("coord","count")
		data[[info[1]]] <<- tab
	}

	else {
		names.4c=NULL
		chroms=NULL
		vps=NULL
		for (i in 1:length(ifiles)) {
			con=file(ifiles[i],"r")
			info=readLines(con,n=1)
			close(con)
			info=unlist(strsplit(info,"\\t"))
			chroms=c(chroms,info[2])
			vps=c(vps,as.numeric(info[3]))
			names.4c=c(names.4c,info[1])
			tab=read.table(ifiles[i],header=FALSE,skip=1)
			colnames(tab)=c("coord","count")
			data[[i]] <<- tab
		}
		chroms=unique(chroms)
		vps=unique(vps)
		if (length(chroms)!=1 | length(vps)!=1) {
			cat("Can only load 4C files with same viewpoint\n")
			data <<- NULL
			return(NULL)
		} else {
			cis.chrom <<- chroms[1]
			vp.pos <<- vps[1]
			names(data) <<- names.4c
		}
	}
	set.conditions(data)
}


#############################################################################
#Save screenshot to file
save.figure=function() {
	fn=tclvalue(tkgetSaveFile(filetypes = "{{Browser screenshots} {.eps}} {{All files} *}"))
	if(nchar(fn)==0) {
		return()
	}
	cat(sprintf("saving figure to: %s\n",fn))
	setEPS()
	postscript(fn)
	min.plot = as.numeric(tclvalue(tkget(min.plot)))
	max.plot = as.numeric(tclvalue(tkget(max.plot)))
	win.plot = as.numeric(tclvalue(tkget(win.plot)))
	win.smooth = as.numeric(tclvalue(tkget(win.smooth)))
	baitname = tclvalue(tkget(baitname))
	plot.ymax = as.numeric(tclvalue(tkget(plot.ymax)))
	plot.4c(data=data,min.plot,max.plot,win.plot,win.smooth,baitname,plot.ymax,settings=settings)
	dev.off()
}

save.bedgraph=function() {
	min.plot = as.numeric(tclvalue(tkget(min.plot)))
	max.plot = as.numeric(tclvalue(tkget(max.plot)))
	win.plot = as.numeric(tclvalue(tkget(win.plot)))
	win.smooth = as.numeric(tclvalue(tkget(win.smooth)))
	
	table=NULL
	if(length(data)>1) {
		for (i in 1:length(settings)) {
			sel = settings[[i]][["indices"]]
			for (j in 1:length(sel)) {
				tmp=data[[sel[j]]]
				tmp=tmp[tmp$coord>vp.pos-1500000 & tmp$coord<vp.pos+1500000,]
				if(is.null(table)) {
					table=tmp
					colnames(table)[2]=names(data)[sel[j]]
				} else {
					colnames(tmp)[2]=names(data)[sel[j]]
					table=merge(table,tmp)
				}
			}
		}
		n=normalizeBetweenArrays(table[,-1])
		for (i in 2:dim(table)[2]) {
			table[,i]=n[,i-1]
		}
	} else {
		table=data[[1]]
		colnames(table)[2]=names(data)[1]
	}

	if(!is.na(min.plot) & !is.na(max.plot)) {
		table=table[table$coord>=min.plot & table$coord<=max.plot,]
	} else if (!is.na(win.plot)) {
		table=table[table$coord>vp.pos-win.plot & table$coord<vp.pos+win.plot,]
	} else {
		cat("Must select either a window size or a full min/max constraint\n")
		return(NULL)
	}

	for (i in 1:length(settings)) {
		tab=NULL
		entries=settings[[i]][["entries"]]
		for (j in 1:length(entries)) {
			tab=cbind(tab,table[,entries[j]])
		}
		bed=data.frame("value"=rowMeans(tab))
		bed$chr=cis.chrom
		bed$start=table$coord
		bed$end = NA
		bed$end[dim(bed)[1]]=bed$start[dim(bed)[1]]+100
		for (j in 1:dim(bed)[1]-1) {
			bed$end[j]=bed$start[j+1]
		}
		bed=bed[,c(2,3,4,1)]
		fn=paste0(settings[[i]]$name,".bedgraph")
		cat(paste0("Saving ",fn,"\n"))
		write.table(bed,fn,col.names=F,row.names=F,quote=F,sep="\t")
	}
}

###########################################################################
#Plot profile

#Gene plot util
parking=function(left,right) {
	y=rep(-1,length(right))
	lengths=right-left
	for (i in order(lengths,decreasing=TRUE)) {
		otherleft=left
		otherleft[i]=NA
		otherright=right
		otherright[i]=NA
		placed=FALSE
		y[i]=0
		while(placed==FALSE) {
			placed=sum((right[i]>otherleft[y==y[i]])&(left[i]<otherright[y==y[i]]),na.rm=TRUE)==0
			if(placed==FALSE) {
				y[i]=y[i]+1
			}
		}
	}
	return(y)
}

#Gene plot sub
plot.genes=function(genome,x.min,x.max) {
	y_plot=parking(genome$Start,genome$End)
	plot(c(x.min,x.max),c(1,-max(y_plot)-0.5),col="white",ylab="",xlab="",fg="white",col.axis="white",xaxs="i",yaxs="i")
	arrowHeads=pretty(x.min:x.max,n=50)
	for(i in 1:dim(genome)[1]) {
		x=c(genome$Start[i],arrowHeads[arrowHeads>genome$Start[i]&arrowHeads<genome$End[i]],genome$End[i])
		if(genome$Strand[i]=="-") {
			arrows(x[2:length(x)],-y_plot[i],x[1:length(x)-1],col="blue",length=0.08)
		} else {
			arrows(x[1:length(x)-1],-y_plot[i],x[2:length(x)],col="blue",length=0.08)
		}
		text(genome$Start[i],-y_plot[i]+0.4,adj=0,labels=genome$Name[i])
	}
}

#Obtain track autoscale levels
get.track.levels=function(epigenome) {
	levels=NULL
	for (i in 1:length(epigenome)) {
		if(epigenome[[i]]$level>0) {
			levels=c(levels,epigenome[[i]]$level)
		}
	}
	levels=unique(levels[order(levels)])
	if(levels[1]!=1 & levels[length(levels)]!=length(levels)) {
		cat("Track levels must be ascending integers with no gaps: 1,2,3,...\n")
		return(NULL)
	}
	plottracklevels=list()
	for (l in 1:length(levels)) {
		for (i in 1:length(epigenome)) {
			if(epigenome[[i]]$level == levels[l]) {
				if(length(plottracklevels)<l) {
					plottracklevels[[l]]=names(epigenome)[i]
				} else {
					plottracklevels[[l]]=c(plottracklevels[[l]],names(epigenome[i]))
				}
			}
		}
	}
	return(plottracklevels)
}

#Track plot sub
plot.tracks=function(chr,x.min,x.max, tracks) {
	tmps=list()
	plotlim=numeric()
	for (i in 1:length(tracks)) {
		tmp=epigenome[[tracks[i]]]$track
		tmp=tmp[as.character(seqnames(tmp))==chr & start(tmp)>=x.min & end(tmp)<=x.max,]
		tmps[[i]]=tmp
		plotlim=max(c(plotlim,score(tmp)))
	}
	for (i in 1:length(tmps)) {
		plot.new()
		plot.window(xlim=c(x.min,x.max),ylim=c(0,plotlim),xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab="",ylab="")
		segments(start(tmps[[i]]),0,end(tmps[[i]]),score(tmps[[i]]),col=epigenome[[tracks[i]]]$color)
		title(main=tracks[i],col.main=epigenome[[tracks[i]]]$color,adj=0,line=-2)
	} 
}

#Plot profile
plot.4c = function(data,min.plot,max.plot,win.plot,win.smooth,baitname,plot.ymax,settings) {
	#Normalize across all used interaction sets if more than one, restrict to 1.5 Mb window
	table=NULL
	if(length(data)>1) {
		for (i in 1:length(settings)) {
			sel = settings[[i]][["indices"]]
			for (j in 1:length(sel)) {
				tmp=data[[sel[j]]]
				tmp=tmp[tmp$coord>vp.pos-1500000 & tmp$coord<vp.pos+1500000,]
				if(is.null(table)) {
					table=tmp
					colnames(table)[2]=names(data)[sel[j]]
				} else {
					colnames(tmp)[2]=names(data)[sel[j]]
					table=merge(table,tmp)
				}
			}
		}
		n=normalizeBetweenArrays(table[,-1])
		for (i in 2:dim(table)[2]) {
			table[,i]=n[,i-1]
		}
	} else {
		table=data[[1]]
		colnames(table)[2]=names(data)[1]
	}

	#Define plot window (min/max overrides fixed window)
	if(!is.na(min.plot) & !is.na(max.plot)) {
		table=table[table$coord>=min.plot & table$coord<=max.plot,]
	} else if (!is.na(win.plot)) {
		table=table[table$coord>vp.pos-win.plot & table$coord<vp.pos+win.plot,]
	} else {
		cat("Must select either a window size or a full min/max constraint\n")
		return(NULL)
	}

	n.plots=1
	if(length(epigenome)>0) {
		for (i in 1:length(epigenome)) {
			if(epigenome[[i]]$level>0) {
				n.plots=n.plots+1
			}
		}
	}
	if(!is.null(genome)) {
		n.plots=n.plots+1
	}
	if(n.plots>1) {
		layout(matrix(1:n.plots,ncol=1,nrow=n.plots),heights=c(20,rep(10,n.plots-1)))
	}
	#first plot - 4C
	x.min=min(table$coord)
	x.max=max(table$coord)
	avg=NULL
	if(length(data)>1) {
		tab=NULL
		entries=settings[[1]][["entries"]]
		for (i in 1:length(entries)) {
			tab=cbind(tab,table[,entries[i]])
		}
		avg=rowMeans(tab)
	} else {
		avg=table[,2]
	}
	avg=caTools::runmean(avg,k=win.smooth,endrule="mean",align="center")
	if (plot.ymax==0) {
		plot.ymax=max(avg)
	}
	main.text = paste0(baitname,"- ",cis.chrom,": ",x.min,"-",x.max,"\n")
	plot(x=table$coord, y=avg,type="l",col=settings[[1]][["color"]],ylab="Running mean QN signal",xlab="Genomic coordinates",xaxs="i",yaxs="i",ylim=c(0,plot.ymax),lwd=2,main=main.text)
	if(length(settings)>1) {
		for (i in 2:length(settings)) {
			tab=NULL
			entries=settings[[i]][["entries"]]
			for (j in 1:length(entries)) {
				tab=cbind(tab,table[,entries[j]])
			}
			avg=rowMeans(tab)
			lines(x=table$coord, y=caTools::runmean(avg,k=win.smooth,endrule="mean",align="center"),col=settings[[i]]$color,lwd=2)
		}
	}
	plotnames=c()
	plotcolors=c()
	for (i in 1:length(settings)) {
		plotnames=c(plotnames,settings[[i]]$name)
		plotcolors=c(plotcolors,settings[[i]]$color)
	}
	legend("topleft",inset=0.05,plotnames,col=plotcolors,lty=1.5,title="4C profile")
	abline(v=vp.pos,col="black")

	#Overlay interactions
	if(length(ints)>0) {
		keep=list()
		for (i in 1:length(ints)) {
			if(ints[[i]]$plot == "1") {
				keep[[ints[[i]]$name]]=ints[[i]]
			}
		}
		if(length(keep)>0) {
			for (i in 1:length(keep)) {
				inttab=keep[[i]]$bed
				inttab=inttab[inttab$chr==cis.chrom & inttab$start>=x.min & inttab$end<=x.max,]
				intcol=keep[[i]]$color
				rect(inttab$start,rep(0,dim(inttab)[1]),inttab$end,rep(plot.ymax,dim(inttab)[1]),border=intcol)
			}
		}
		intcols=c()
		for (i in 1:length(keep)) {
			intcols=c(intcols,keep[[i]]$color)
		}
		legend("left",inset=0.05,names(keep),col=intcols,lty=1.5,title="Interactions")
	}
	n.plots=n.plots-1
	#second plot - gene track
	if(!is.null(genome)) {
		genome = genome[genome$Chr==cis.chrom & genome$Start>=x.min & genome$End<=x.max,]
		if(dim(genome)[1]>0) {
			plot.genes(genome,x.min,x.max)
		}
	n.plots=n.plots-1
	}
	#other plots - epigenomic tracks, with level-specific autoscaling
	if(n.plots>0) {
		plottracklevels <<- get.track.levels(epigenome)
		for (l in 1:length(plottracklevels)) {
			plot.tracks(cis.chrom,x.min,x.max,plottracklevels[[l]])
		}
	}
}

#Wrapper to set up plotting function and output to gui
refresh.gui = function() {
	if(is.null(data)) {
		cat("Need to open 4C file(s) first - see menu file\n")
		return(NULL)
	}
	min.plot = as.numeric(tclvalue(tkget(min.plot)))
	max.plot = as.numeric(tclvalue(tkget(max.plot)))
	win.plot = as.numeric(tclvalue(tkget(win.plot)))
	win.smooth = as.numeric(tclvalue(tkget(win.smooth)))
	baitname = tclvalue(tkget(baitname))
	plot.ymax = as.numeric(tclvalue(tkget(plot.ymax)))
	hscale = 1.99
	vscale = 1.99
	if(is.null(browse$env$plot)) {
		options(show.error.messages = FALSE) #
		browse$env$plot = try(tkrplot(browse, plot.4c(data=data,min.plot,max.plot,win.plot,win.smooth,baitname,plot.ymax,settings=settings),hscale=hscale,vscale=vscale),silent=TRUE)
		tkgrid(browse$env$plot)
	}
	else {
		options(show.error.messages = FALSE)
		try(tkrreplot(browse$env$plot, plot.4c(data=data,min.plot,max.plot,win.plot,win.smooth,baitname,plot.ymax,settings=settings),hscale=hscale,vscale=vscale),silent=TRUE)
		tkgrid(browse$env$plot)
	}
}

########################################################################################################################

#GUI

main = tktoplevel()
browse = tktoplevel()
tktitle(main) = "Main"
tktitle(browse) = "Browser"
main$env$menu = tk2menu(main)
tkconfigure(main,menu=main$env$menu)

#Main menu (File -> New,Save Image,Quit)
main$env$menuFile = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="File",menu=main$env$menuFile)
tkadd(main$env$menuFile, "command", label = "New", command=load.4c)
tkadd(main$env$menuFile, "command", label = "Save Image", command=save.figure)
tkadd(main$env$menuFile, "command", label = "Save Bedgraph(s)", command=save.bedgraph)
tkadd(main$env$menuFile, "command", label = "Quit", command=function() tkdestroy(main))

#Conditions menu (-> Set Conditions,Plot Conditions)
main$env$menuCon = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="Conditions", menu=main$env$menuCon)
tkadd(main$env$menuCon, "command", label="Set Conditions", command=function() set.conditions(data))
tkadd(main$env$menuCon, "command", label="Plot Conditions", command=manage.set)

#Interactions menu (-> Load Interactions,Manage Interaction Plots)
main$env$menuInt = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="Interactions", menu=main$env$menuInt)
tkadd(main$env$menuInt, "command", label="Load Interactions", command=load.ints)
tkadd(main$env$menuInt, "command", label="Manage Interaction Plots", command=manage.ints)

#Tracks menu (-> Load Genes, Load Tracks, Manage Tracks)
main$env$menuTrk = tk2menu(main$env$menu, tearoff=FALSE)
tkadd(main$env$menu, "cascade", label="Tracks", menu=main$env$menuTrk)
tkadd(main$env$menuTrk, "command", label="Load Genes", command=load.genes)
tkadd(main$env$menuTrk, "command", label="Load Tracks", command=load.epigenome)
tkadd(main$env$menuTrk, "command", label="Manage Tracks", command=manage.tracks)


#Plot window and bait choice
labmin = tklabel(main,text="start coordinate")
min.plot <<- tkentry(main,width=20, textvariable=tclVar("NA"))
tkgrid(labmin,row=1,column=1,sticky="e")
tkgrid(min.plot,row=1,column=2,sticky="w")

labmax = tklabel(main,text="end coordinate")
max.plot <<- tkentry(main,width=20,textvariable=tclVar("NA"))
tkgrid(labmax,row=2,column=1,sticky="e")
tkgrid(max.plot,row=2,column=2,sticky="w")

labwin = tklabel(main,text="plot window")
win.plot <<- tkentry(main, width=10, textvariable=tclVar("300000"))
tkgrid(labwin, row=3, column=1, sticky="e")
tkgrid(win.plot, row=3, column=2, sticky="w")

labsmooth = tklabel(main,text="smooth window")
win.smooth <<- tkentry(main, width=10, textvariable=tclVar("21"))
tkgrid(labsmooth, row=4, column=1, sticky="e")
tkgrid(win.smooth, row=4, column=2, sticky="w")

labbait = tklabel(main,text="bait name")
baitname <<- tkentry(main, width=10, textvariable=tclVar("Bait"))
tkgrid(labbait, row=5, column=1, sticky="e")
tkgrid(baitname, row=5, column=2, sticky="w")

labymax = tklabel(main,text="max y plot")
plot.ymax <<- tkentry(main,width=10,textvariable=tclVar("0"))
tkgrid(labymax, row=6, column=1, sticky="e")
tkgrid(plot.ymax, row=6, column=2, sticky="w")

ok=tkbutton(main,text="OK",command=refresh.gui)
tkgrid(ok,row=7,column=1)
