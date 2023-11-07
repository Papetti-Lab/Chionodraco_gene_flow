# ------------------------------------------ #
# Modified functions to plot TreeMix results #
# ------------------------------------------ #

## Additional arguments in the function 'custom_plot_tree' to make plots more customisable and clearer visualisation
# lwd.arrow =      line width for the arrow
# cex.x.label =    font size for the x axis label
# cex.x.values =   font size for the x axis values
# line.x =         distance of the x label from the axis
# padj.x =         distance of the x values from the axis
# cex.colorbar =   font size for colour bar labelling
# cex.mse =        font size for s.e. bar labelling
# se.y =           value to adjust the position of the se bar and labelling on the y axis
# cb.y =           value to adjust the position of the colour bar and labelling on the y axis


custom_plot_tree = function(stem, o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1,
lwd.arrow=1, cex.x.label=1, cex.x.values=1, line.x=1, padj.x=0.5, cex.colorbar=1, cex.mse=1, se.y=0, cb.y=0){
  d = paste(stem, ".vertices.gz", sep = "")
  e = paste(stem, ".edges.gz", sep = "")
  se = paste(stem, ".covse.gz", sep = "")
  d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
  e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
  if (!is.na(o)){
    o = read.table(o, as.is = T, comment.char = "", quote = "")
  }
  e[,3] = e[,3]*e[,4]
  e[,3] = e[,3]*e[,4]
  
  se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
  m1 = apply(se, 1, mean)
  m = mean(m1)
  #m = 0
  for(i in 1:length(flip)){
    d = flip_node(d, flip[i])
  }
  d$x = "NA"
  d$y = "NA"
  d$ymin = "NA"
  d$ymax = "NA"
  d$x = as.numeric(d$x)
  d$y = as.numeric(d$y)
  d$ymin = as.numeric(d$ymin)
  d$ymax = as.numeric(d$ymax)
  
  d = set_y_coords(d)
  d = set_x_coords(d, e)
  print(d)
  d = set_mig_coords(d, e)
  custom_plot_tree_internal(d, e, o = o, cex = cex, xmin = xmin, disp = disp, plus = plus, arrow = arrow, ybar = ybar, mbar = mbar, mse = m, scale = scale, plotmig = plotmig, plotnames = plotnames, lwd = lwd, font = font,
lwd.arrow=lwd.arrow, cex.x.label=cex.x.label, cex.x.values=cex.x.values, line.x=line.x, padj.x=padj.x, cex.colorbar=cex.colorbar, cex.mse=cex.mse, se.y=se.y, cb.y=cb.y)
  return(list( d= d, e = e))
}

custom_plot_tree_internal = function(d, e, o = NA, cex = 1, disp = 0.005, plus = 0.005, arrow = 0.05, ybar = 0.01, scale = T, mbar = T, mse = 0.01, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1,
lwd.arrow=1, cex.x.label=1, cex.x.values=1, line.x=1, padj.x=0.5, cex.colorbar=1, cex.mse=1, se.y=se.y, cb.y=cb.y){
  plot(d$x, d$y, axes = F, ylab = "", xlab = "", xlim = c(xmin, max(d$x)+plus), pch = "")
  title(xlab="Drift parameter", line=line.x, cex.lab=cex.x.label)
  axis(1, cex.axis = cex.x.values, padj = padj.x)
  mw = max(e[e[,5]=="MIG",4])
  mcols = rev(heat.colors(150))
  for(i in 1:nrow(e)){
    col = "black"
    if (e[i,5] == "MIG"){
      w = floor(e[i,4]*200)+50
      if (mw > 0.5){
        w = floor(e[i,4]*100)+50
      }
      col = mcols[w]
      if (is.na(col)){
        col = "blue"
      }
    }
    v1 = d[d[,1] == e[i,1],]
    v2 = d[d[,1] == e[i,2],]
    if (e[i,5] == "MIG"){
      if (plotmig){
        arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, col = col, length = arrow, lwd = lwd.arrow)
      }
    }
    else{
      lines( c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y), col = col, lwd = lwd)
    }
  }
  tmp = d[d[,5] == "TIP",]
  print(tmp$x)
  print(disp)
  if ( !is.na(o)){
    for(i in 1:nrow(tmp)){
      tcol = o[o[,1] == tmp[i,2],2]
      if(plotnames){
        #print(tmp[i,2])
        text(tmp[i,]$x+disp, tmp[i,]$y, labels = tmp[i,2], adj = 0, cex = cex, col  = tcol, font = font, xpd = T)
      }
    }
  }
  else{
    if (plotnames){
      text(tmp$x+disp, tmp$y, labels = tmp[,2], adj = 0, cex = cex, font = font, xpd = T)
    }
  }
  if (scale){
    print (paste("mse", mse))
    lines(c(0, mse*10), c(ybar + se.y, ybar + se.y))
    text( 0, ybar + 0.03 + se.y, lab = "10 s.e.", adj = 0, cex  = cex.mse)
    lines( c(0, 0), c( ybar - 0.01 + se.y, ybar+0.01 + se.y))
    lines( c(mse*10, mse*10), c(ybar- 0.01 + se.y, ybar+ 0.01 + se.y))
  }
  if (mbar){
    mcols = rev( heat.colors(150) )
    mcols = mcols[50:length(mcols)]
    ymi = ybar+0.15
    yma = ybar+0.35
    l = 0.2
    w = l/100
    xma = max(d$x/20)
    rect( rep(0, 100), (ymi+(0:99)*w) + cb.y, rep(xma, 100), (ymi+(1:100)*w) + cb.y, col = mcols, border = mcols)
    text(xma+disp, ymi + cb.y, lab = "0", adj = 0, cex = cex.colorbar)
    if ( mw >0.5){ text(xma+disp, yma + cb.y, lab = "1", adj = 0, cex = cex.colorbar)}
    else{
      text(xma+disp, yma + cb.y, lab = "0.5", adj = 0, cex = cex.colorbar)
    }
    text(0, yma+0.09 + cb.y, lab = "Migration", adj = 0 , cex = cex.colorbar)
    text(0, yma+0.05 + cb.y, lab = "weight", adj = 0 , cex = cex.colorbar)
  }	
}


################################################################################


## Additional arguments in the function 'custom_plot_tree' to make plots more customisable and clearer visualisation
# font      font for labelling
# cex.bar   font size for colour bar labelling
# y.pos     position of the labelling on the y direction


custom_plot_resid = function(stem, pop_order, min = -0.009, max = 0.009, cex = 1, usemax = T, wcols = "r",
                         font=1, cex.bar=1, y.pos=1){
  c = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  m = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  names(c) = rownames(c)
  names(m) = rownames(m)
  o = read.table(pop_order, as.is = T, comment.char = "", quote = "")
  se = read.table(gzfile(paste(stem, ".covse.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  mse = apply(se, 1, mean)
  mse = mean(mse)
  print(mse)	
  c = c[order(names(c)), order(names(c))]
  m = m[order(names(m)), order(names(m))]
  tmp = c -m 
  #tmp = m - c
  #tmp = (m-c)/m
  #print(tmp)
  toplot = data.frame(matrix(nrow = nrow(tmp), ncol = ncol(tmp)))
  for(i in 1:nrow(o)){
    for( j in 1:nrow(o)){
      #print(paste(o[i,1], o[j,1]))
      if (o[i,1] %in% names(tmp) ==F){
        print(paste("not found", o[i,1]))
      }
      if (o[j,1] %in% names(tmp) ==F){
        print(paste("not found", o[j,1]))
      }
      toplot[i, j] = tmp[which(names(tmp)==o[i,1]), which(names(tmp)==o[j,1])]
    }
  }
  #print(toplot)
  if (usemax){
    m1 = max(abs(toplot), na.rm = T)
    max = m1*1.02
    min = -(m1*1.02)	
  }
  print("here")
  names(toplot) = o[,1]
  toreturn = custom_plot_resid_internal(toplot, max = max, min = min, wcols = wcols, mse = mse, o = o, cex = cex,
                                    font=font, cex.bar=cex.bar, y.pos=y.pos)
  return(toreturn)
}


custom_plot_resid_internal = function(d, o = NA, max = 0.009, min = -0.009, cex =0.5, wcols = "rb", mse = NA,
                                  font=1, cex.bar=cex.bar, y.pos=y.pos){
  npop = nrow(d)
  width = 1/npop
  height = 1/npop
  colors = brewer.pal(9, "Spectral")
  colors = c("red", "orange","yellow", "white", "green", "blue", "black")
  pal = colorRampPalette(colors)
  ncol = 80
  cols = pal(ncol)
  plot("NA", xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = "", ylab = "")
  for (i in 1:npop){
    for( j in 1:i){
      v = d[i,j]
      print(paste(i, j, v))
      col= "white"
      if (v < 0){
        if (wcols == "rb"){
          col = rgb(0, 0, 1, v/min)
        }
        else{
          #col = rgb(0, 0, 1, 0.1+0.9*(v/min))
          col = cols[ncol/2-floor( (v/min)*(ncol/2))]
          #col = "white"
        }
      }
      else{
        if (wcols == "rb"){
          col = rgb(1, 0, 0, v/max)
        }
        else{
          #col = rgb(1, 0, 0, 0.1+0.9*(v/max))
          col = cols[ncol/2+ceiling((v/max)*(ncol/2))]
        }
      }
      xmin = j/npop - 1/npop
      xmax = j/npop
      ymin = 1-(i/npop)
      ymax = 1-(i/npop)+1/npop
      rect(xmin, ymin, xmax, ymax, col = col, border = "black")
    }
    tcol = "black"
    tmp = o[o[,1] == names(d)[i],]
    if (length(tmp) != 1){
      tcol = tmp[1,2]
    }
    mtext(names(d)[i], side = 2, at = 1-i/npop+0.5/npop, las = 1, cex = cex, col = tcol, font=font, xpd=T)
    #mtext(names(d)[i], side = 1, at =  i/npop-0.5/npop, las = 3, cex = cex, col = tcol, font=font)
    text(i/npop - 0.5/npop, y.pos, labels = names(d)[i], adj = 1, srt = 40, cex = cex, col = tcol, font=font, xpd=T)
    
  }
  if ( !is.na(mse)){
    ymi = 0.5
    yma = 0.9
    w = (yma-ymi)/ncol
    xma = 0.80
    lmi = round(min/mse, digits = 1)
    lma = round(max/mse, digits = 1)
    print(cols)
    print(ymi+(0:ncol)*w)
    rect( rep(0.75, ncol), ymi+(0:(ncol-1))*w, rep(xma, ncol), ymi+(1:ncol)*w, col = cols, border = cols)
    text(xma+0.01, ymi, lab = paste(lmi, "SE"),  adj = 0, cex = cex.bar)
    text(xma+0.01, yma, lab = paste(lma, "SE"), adj = 0, cex = cex.bar)
    
  }
  return(d)
  #image(as.matrix(d), col = cols)
}


