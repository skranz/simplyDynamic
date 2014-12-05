
showcontours <- function(fun,xlims, ylims, resol=50,add=F, colors=c('red', 'blue')) {
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1, add=add, col=colors[1]);
  contour(x[,1],y[1,],z2, add=T, col=colors[2]); 
}

nullclines <- function(fun,xlims, ylims, resol=100, add=F,colors=c('green', 'blue'),...) {
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1,levels=c(0), add=add, col=colors[1],...);
  contour(x[,1],y[1,],z2,levels=c(0), add=T, col=colors[2],...); 
}

phasearrows <- function(fun,xlims,ylims,resol=10, col='grey', add=F) {
  if (add==F) {
    plot(1,xlim=xlims, ylim=ylims, type='n');
  }
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=T, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=F, resol, resol);
  npts <- resol*resol;
  xspace <- abs(diff(xlims))/(resol*5);
  yspace <- abs(diff(ylims))/(resol*5);
  x <- x + matrix(runif(npts, -xspace, xspace),resol,resol);
  y <- y + matrix(runif(npts, -yspace, yspace),resol,resol);
  z <- fun(x,y);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  maxx <- max(abs(z1));
  maxy <- max(abs(z2));
  dt <- min( abs(diff(xlims))/maxx, abs(diff(ylims))/maxy)/resol;
  lens <- sqrt(z1^2 + z2^2);
  lens2 <- lens/max(lens); 

  arrows(c(x), c(y), c(x+dt*z1/((lens2)+.1)), c(y+dt*z2/((lens2)+.1)),length=.04, col=col);
}

get.lhs.and.rhs = function(txt) {
  txt = sep.lines(txt, collapse=";")
  txt = sep.lines(txt, collapse="\n")

  mat = str_trim(do.call(rbind,str.split(txt,"=")))
  colnames(mat) = c("lhs","rhs")  
  mat
}

process.model.text = function(txt) {
  library(skUtils)
  txt = sep.lines(txt, collapse=";")
  txt = sep.lines(txt, collapse="\n")
  txt = str_trim(txt)
  remove = txt=="" | substring(txt,1,1)=="#"
  txt = txt[!remove]
  
  # Extract equations and variables
  mat = str_trim(do.call(rbind,str.split(txt,"=")))
  names = mat[,1]

  eq.rows = substring(names,1,4) == "dot." 
  var = str.remove.ends(mat[eq.rows,1],4)
  eq = mat[eq.rows,2]
  var.list = paste(var,collapse=",")
  
  # Extract parameter
  mat = mat[!eq.rows,]
  rownames(mat)=mat[,1]
  param.names = setdiff(mat[,1],var)
  param = as.numeric(mat[param.names,2])
  names(param) = param.names
  init = as.numeric(mat[var,2])
  names(init) = var
  
  # Evaluate param in an environment
  env = new.env(.GlobalEnv)
  param.run = paste(param.names, " <- ", param, ";")
  eval(parse(text=param.run),env)

  # Make function for plotting arrows
  code = paste(
    'function(',var.list,') {
      c(',paste0('"',var,'"=',eq,collapse=","),')
    }'
  )
  fun = eval(parse(text=code))
  environment(fun) = env
  
  # Make model for solving ODE

  code = paste0('function (time, y, parms) {
    with(as.list(c(y, parms)), {
      list(c(',paste0("d",var, " = ", eq, collapse=","),'))
    })
  }')
  model = eval(parse(text=code))
  environment(model) = .GlobalEnv

  code = paste0('function (time, y, parms) {
    with(as.list(c(y, parms)), {
      list(c(',paste0("d",var, " = ", eq, " + ", var, collapse=","),'))
})
}')
  new.val.model = eval(parse(text=code))
  environment(new.val.model) = .GlobalEnv
  
  
  return(list(n=length(var),var=var,param=param, init=init, model=model, new.val.model=new.val.model, fun=fun))              
}

examples.process.model.text = function() {
  txt = "
# Hi
  dot.X = pmin(-X,X*(a-b*Y))
  dot.Y = pmin(-Y,-Y*(c-d*X))
# Var
  X = 1;Y=2
  a = 1
  b = 0.1
  c=0.1
  d=0.1
"
  mod = process.model.text(txt)
  
  out <- ode(mod$init, seq(10000,11000,by=1), mod$model, mod$param)
  if (mod$n == 2) {
    var = mod$var
    X = out[,2]
    Y = out[,3]
    xlim = range(X)
    ylim = range(Y)
    plot(X,Y,xlim=xlim,ylim=ylim,xlab=var[1],ylab=var[2],type="l",col="red",lwd=2)
    
    phasearrows(mod$fun,xlim,ylim,add=TRUE)
    nullclines(mod$fun,xlim,ylim,add=TRUE,lwd=2)
    
  } 
}

examples.phasearrows = function() {
  a = 1; b = 0.1; c=0.1; d=0.1
  fun = function(X,Y) {
    dX <- X*(a-b*Y)
    dY <- -Y*(c-d*X)
    c(dX,dY)
  }
  
  
  phasearrows(fun,c(-10,100),c(-10,100),20, col = "grey")
  par(mfrow=c(1,1))
}

plot.speed.points = function(x,y,pch=19,rows=1:length(x),cex=1.5,nh=7,...) {
  x = x[rows]; y = y[rows]
  n = length(x)
  head.col = hsv(1/6, c(1,seq(0.9,0.5, length=nh-1)),c(1,seq(1,0.8, length=nh-1)))
  col = c(rep("grey",max(0,n-nh)),rev(head.col[1:min(n,nh)]))
  for (i in length(x):1)
    points(x[i],y[i],col=col[i],pch=pch,cex=cex)
  #points(x,y,col=col,pch=pch,...)
}

animate.speed.points = function(mod,out, sleep=0.1, nh = 20) {
  # Open separate plot window
  
  win.graph()
  var = mod$var
  X = out[,2]
  Y = out[,3]
  xlim = range(X)
  ylim = range(Y)
  plot(X,Y,xlim=xlim,ylim=ylim,xlab=var[1],ylab=var[2],type="l",col="grey",lwd=2)
  phasearrows(mod$fun,xlim,ylim,add=TRUE)
  nullclines(mod$fun,xlim,ylim,add=TRUE,lwd=2)
   
  for (i in 1:length(X)) {
    rows = max(1,i-nh-1):i
    plot.speed.points(X,Y,rows=rows,nh=nh)
    Sys.sleep(sleep)
  }
  
  
}

examples.plot.speed.points = function() {
  txt = "
# Hi
# dot.X = pmax(-X,(X*(a-b*Y) )*0.1+rnorm(length(X),0.1,0.1))
# dot.Y = pmax(-Y,((-Y*(c-d*X))*0.1-Y*0.5+rnorm(length(X),0,0.2)))
  
  dot.X = 0.1+(X*(a-b*Y))
  dot.Y = (Y*(d*X-c))

  # Var
  X = 5;Y=2
  a = 1
  b = 0.1
  c=0.15
  d=0.1
  "
  mod = process.model.text(txt)
  out <- ode(mod$init, seq(0,100,by=0.05), mod$model, mod$param)
  
  #out <- ode(mod$init, seq(0,10000,by=1), mod$new.val.model, mod$param, method="iteration")
  animate.speed.points(mod,out,sleep=.05)
  var = mod$var
  X = out[,2]
  Y = out[,3]
  xlim = range(X)
  ylim = range(Y)
  plot(X,Y,xlim=xlim,ylim=ylim,xlab=var[1],ylab=var[2],type="l",col="red",lwd=2)
  
  rows = 1:10
  x = X
  y = Y
  plot.speed.points(x,y,rows=1:10)
  
  phasearrows(mod$fun,xlim,ylim,add=TRUE)
  nullclines(mod$fun,xlim,ylim,add=TRUE,lwd=2)
 
  
}