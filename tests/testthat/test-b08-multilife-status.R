test_that("C17 rejects incompatible complete lifetimes before aggregation", {
  for(status in c("joint","last_survivor")) for(tol in c(1e-10,1e-3)) {
    for(time in c(2-1e-12,3,3+1e-12)) {
      d<-data.frame(sim_id=c(1,1),life_id=c("x","y"),Kx=c(2,5),Tx=c(time,5.5))
      expect_error(mc_multilife_status(d,status=status,tol=tol),"K <= T < K")
    }
  }
})

test_that("valid boundary lifetimes retain exact min and max status books", {
  d<-data.frame(sim_id=c(3,3,1,1,2,2),life_id=rep(c("x","y"),3),
                Kx=c(2,5,1,4,3,2),Tx=c(2,5.5,2-1e-12,4.25,3.5,2.75))
  for(status in c("joint","joint_life","last","last_survivor","last_survivor_life")) {
    fun<-if(status %in% c("joint","joint_life")) min else max
    rows<-list(1:2,3:4,5:6)
    kk<-vapply(rows,function(j) fun(d$Kx[j]),0)
    tt<-vapply(rows,function(j) fun(d$Tx[j]),0)
    z<-mc_multilife_status(d,status=status)
    expect_identical(z$sim_id,c(3,1,2))
    expect_equal(z$K_status,kk)
    expect_equal(z$T_status,tt)
    expect_true(all(z$K_status==floor(z$T_status)))
    expect_equal(mc_insurance(z,i=0,benefit=10,col_K="K_status",col_T="T_status")$pv_benefit,
                 rep(10,3))
    expect_equal(mc_annuity(z,i=0,payment=1,timing="due",col_K="K_status",col_T="T_status")$pv_annuity,
                 kk+1)
  }
  # tol still checks integer K, without replacing a valid complete lifetime.
  d$Kx[1]<-2+1e-12
  z<-mc_multilife_status(d,tol=1e-10)
  expect_equal(z$K_status,c(2,1,2))
  expect_equal(z$T_status,c(2,2-1e-12,2.75))
})

test_that("large status lifetimes do not overflow integer output storage", {
  d<-data.frame(sim_id=c(1,1),life_id=c("x","y"),
                Kx=c(3e9,3e9+2),Tx=c(3e9+.5,3e9+2.75))
  for(status in c("joint","last_survivor")) {
    fun<-if(status=="joint") min else max
    out<-mc_multilife_status(d,status=status)
    expect_equal(out$K_status,fun(d$Kx))
    expect_equal(out$T_status,fun(d$Tx))
    expect_false(anyNA(out$K_status))
  }
  d$Kx<-c(2,3); d$Tx<-c(2.5,3.5)
  expect_type(mc_multilife_status(d)$K_status,"integer")
})

test_that("simulation identifiers cannot collide with output fields", {
  for(nm in c("K_status","T_status","n_lives","status")) {
    d<-data.frame(id=c(1,1),life_id=c("x","y"),Kx=c(2,5),Tx=c(2.5,5.5))
    names(d)[1]<-nm
    expect_error(mc_multilife_status(d,col_sim=nm),"reserved|output")
  }
  d<-data.frame(id=c("a","a","b","b"),life_id=rep(c("x","y"),2),
                Kx=c(1,3,2,4),Tx=c(1.2,3.5,2.4,4.6))
  out<-mc_multilife_status(d,col_sim="id")
  expect_identical(out$id,c("a","b"))
  expect_identical(out$sim_id,out$id)
  expect_identical(out$sim,out$id)
})
