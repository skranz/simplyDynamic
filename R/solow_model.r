
examples.solow.model = function() {
  m = dynamic.model({
    s = 0.1
    L = 1
    delta = 0.1
    
    K = 2
    
    K[t] = I[t-1] +(1- delta) * K[t-1]
    Y[t] = sqrt(K[t]) * sqrt(L)
    I[t] = s*Y[t]
    
  })

  # Two groups of population with unequal saving rates
  m = dynamic.model({
    # We have two groups of the population with equal size
    s1 = 0.01
    s2 = 0.2
    K.share = 0.8
    L.share = 1-K.share
    share1 = 0.9
    share2 = 1-share1
    tax = 0.2
    
    L = 1
    delta = 0.1
    
    K1 = 1
    K2 = 0.1
    
    
    K1[t] = I1[t-1] +(1- delta) * K1[t-1]
    K2[t] = I2[t-1] +(1- delta) * K2[t-1]
    K[t] = K1[t]+K2[t]
    
    Y[t] = (K[t]^K.share) * (L^(1-K.share))
    Y1.gross[t] = Y[t] *L.share*share1 + Y[t] *K.share*K1[t]/K[t]
    Y2.gross[t] = Y[t] *L.share*share2 + Y[t] *K.share*K2[t]/K[t]
    
    Y1[t] = Y1.gross[t]*(1-tax) + Y[t]*tax*share1
    Y2[t] = Y2.gross[t]*(1-tax) + Y[t]*tax*share2
    
    
    I1[t] = s1*Y1[t]
    I2[t] = s2*Y2[t]
    
    c1[t] = (Y1[t]-I1[t])/share1
    c2[t] = (Y2[t]-I2[t])/share2
    
    
  })

  sim = simulate.model(m,T=500,par = list(K.share=0.5, share1=0.99,tax=0.5))
  plot(sim$t, sim$Y2)
  lines(sim$t,sim$Y1, col="blue")
  sim$K2 / (sim$K1+sim$K2)
  sim$Y2 / sim$Y1
  sim$c2 / sim$c1

  
  sim = simulate.model(m, par=list(K=5))
  plot(sim$t, sim$Y)
  
  solve.steady.state(m)
  
  plot(sim$K,sim$Y)
}
