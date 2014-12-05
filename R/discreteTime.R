dynamic.model = function(code, parent.env = baseenv()) {
  code = substitute(code)
  
  restore.point("dynamic.model")
  calls = code[-1]
  co.li = lapply(calls, parse.model.code.line)
  library(data.table)
  df = as.data.frame(rbindlist(co.li))
   
  vars = unique(df$sym[df$type=="eq"])
  pars = setdiff(unique(df$sym), vars)
  syms = c(pars, vars) 
  
  init.env = new.env(parent=parent.env)
  init.calls = as.list(calls[df$type=="par"])
  names(init.calls) = df$sym[df$type=="par"]
  for (call in init.calls) {
    eval(call,init.env)    
  }
  
  
  rows = df$type=="par" & df$var %in% vars
  df$type[rows] = "init"
  
  eq.calls = as.list(calls[df$type=="eq"])
  names(eq.calls) = df$sym[df$type=="eq"]
  
  model = list(par=pars, var=vars,symbols=syms, init.calls=init.calls, init.env=init.env, eq.calls=eq.calls, df=df, parent.env=parent.env)
  
  model 
}

parse.model.code.line = function(call) {
  lhs = get.lhs(call)
  vi = extract.var.with.index(lhs, as.character=TRUE)
  
  if (is.null(vi$index)) {
    type = "par"
    vi$index = ""
  } else {
    type = "eq"
  }
  list(type=type, sym=vi$var, index=vi$index)
}

#' Simulate a model
simulate.model = function(m,T=100,name=paste0("sim_",sample.int(1e9,1)), par=NULL) {
  restore.point("simulate.model")

  res = init.sim.env(m=m,T=T, par=par)
  sim.env = res$sim.env; init=res$init
  
  # run simulation for all periods
  for (t in 2:T) {
    sim.env$t = t
    for (eq in m$eq.calls) {
      eval(eq,sim.env)
    }
  }
  
  library(dplyr)
  li = lapply(m$var, get, pos=sim.env)
  names(li) = m$var
  sim = do.call(data_frame,li)
  sim = cbind(data_frame(t=1:T),sim)
  attr(sim, "init") = init
  attr(sim, "m") = m
  attr(sim,"name") = name
  sim
}

init.sim.env = function(m,T, par=NULL) {
  restore.point("init.sim.env")
  # Update init env if par is given
  if (length(par)>0) {
    init.env = new.env(parent=m$parent.env)
    init.calls = m$init.calls
    for (p in names(init.calls)) {
      if (p %in% names(par)) {
        init.env[[p]] = par[[p]]
      } else {
        call = init.calls[[p]]
        eval(call,init.env)
      }
    }
  } else {
    init.env = m$init.env
  }
  
  sim.env = new.env(parent=m$parent.env)
  init = as.list(init.env)
  copy.into.env(source=init.env, dest=sim.env,names = m$par)
  
  var = m$var[[1]]  
  for (var in m$var) {
    if (exists(var, init.env, inherits=FALSE)) {
      start = get(var,init.env)
      vals = vector(class(start)[1],T)
      vals[1] = start
    } else {
      vals = numeric(T)
      vals[1] = NA
    }
    sim.env[[var]] = vals
  }
  
  # Variables that are not initialized will be inialized with their
  # dynamic formula
  sim.env$t = 1
  no.init.var = setdiff(m$var,ls(init.env)) 
  for (eq in m$eq.calls[no.init.var]) {
    eval(eq,sim.env)    
  }

  return(list(sim.env=sim.env,init=init))
}

dyn.to.static.eq = function(eq) {
  
  remove.index = function(call) {
    if (is.name(call)) return(call)
    if (as.character(call[[1]])=="[") {
      return(call[[2]])  
    }
    return(call)
  }
  static.eq = recursively.replace(eq, remove.index)
  static.eq
}

solve.steady.state = function(m) {
  library(nleqslv)

  
  static.eq = lapply(m$eq.calls, dyn.to.static.eq)
  F.li = lapply(static.eq, function(eq) {
    substitute(lhs-(rhs),list(lhs=get.lhs(eq), rhs=get.rhs(eq)))
  })
  F = make.call("c", F.li, use.names=FALSE)
  
  # Substitute parameters for their numerical values
  par.list = as.list(m$init.env)[m$par]
  F = substitute.call(F,par.list)
  
  # Assign variables from values
  assign.li = lapply(seq_along(m$var), function(i) {
    substitute(var<-values[i], list(var=as.name(m$var[i]),i=i))
  })
  
  body = make.call("{", c(assign.li, list(F)))
  fn = function(values) {}
  body(fn) <- body
  fn  
  
  sim.env = init.sim.env(m,T=1)$sim.env
  start = unlist(as.list(sim.env)[m$var])
  
  res = nleqslv(start,fn)
  x = res$x
  names(x) = m$var
  x
}
