test_that("reserve benefits retain original rates at multiple durations", {
  for (force in c(720, -1000)) {
    amount <- if (force > 0) 1e300 else 1e-300
    out <- mc_reserve(data.frame(Kx=0,Tx=.75),t=c(0,.25,.5),i=force,
      i_type="force",P=0,benefit=amount)
    expected <- exp(log(amount)-force*(1-out$t))
    expect_equal(out$Z_t/expected,rep(1,3),tolerance=1e-12)
    expect_identical(out$future_pv_premiums,rep(0,3))
    expect_equal(out$L_t/expected,rep(1,3),tolerance=1e-12)
  }
})

test_that("premium flows are not reconstructed from a lost annuity basis", {
  for (unit in c("annualized","per_payment")) for (force in c(720,-1000)) {
    amount <- if (force > 0) 1e300 else 1e-300
    basis <- 1/amount
    out <- mc_reserve(data.frame(Kx=1,Tx=1.5),t=c(0,.5),i=force,i_type="force",
      P=amount,premium_unit=unit,payment=basis,benefit=0,premium_timing="immediate")
    expected <- exp(log(amount)-force*(1-out$t))
    expect_equal(out$future_pv_premiums/expected,c(1,1),tolerance=1e-12)
    expect_equal(out$L_t/expected,c(-1,-1),tolerance=1e-12)
  }
  out <- mc_reserve(data.frame(Kx=1,Tx=1.5),t=0,i=1000,i_type="force",
    P=1e300,premium_unit="annuity_scale",payment=1e300,benefit=0,
    premium_timing="immediate")
  expect_equal(out$future_pv_premiums/exp(2*log(1e300)-1000),1,tolerance=1e-12)
})

test_that("a zero premium remains zero even when its annuity basis overflows", {
  out <- mc_reserve(data.frame(Kx=1,Tx=1.5),t=c(0,1),i=-1000,i_type="force",
                    P=0,benefit=0)
  expect_identical(out$future_pv_premiums,c(0,0))
  expect_identical(out$L_t,c(0,0))
})

test_that("internal premium estimation uses complete issue pairs", {
  d <- data.frame(Kx=c(1,NA,2),Tx=c(1.5,5.5,2.5))
  out <- mc_reserve(d,t=0,i=0,benefit=100,timing="moment_of_death")
  expect_equal(out$premium_annuity_scale,rep(40,3))
  expect_equal(out$L_t,c(20,NA,-20))
  expect_true(all(out$premium_was_estimated))
})

test_that("C20 makes internal estimation independent of input premium units", {
  d <- data.frame(Kx=c(1,2),Tx=c(1.5,2.5))
  for (unit in c("auto","annualized","per_payment","annuity_scale")) {
    out <- mc_reserve(d,t=0,i=0,benefit=100,k=2,payment=7,premium_unit=unit)
    expect_equal(out$premium_per_payment,c(25,25),tolerance=1e-12)
    expect_equal(out$premium_annualized,c(50,50),tolerance=1e-12)
    expect_lt(abs(mean(out$L_t)),1e-12)
  }
})

test_that("C19 partitions term and survival benefits at exact maturity", {
  d <- data.frame(Tx=c(2-1e-12,2,2+1e-12,3,NA))
  d$Kx <- floor(d$Tx)
  parts <- lapply(c("term","pure_endowment","endowment"),function(type)
    mc_reserve(d,t=c(0,1),i=.04,P=0,benefit=100,type=type,n=2,
               timing="moment_of_death"))
  expect_equal(parts[[1]]$Z_t+parts[[2]]$Z_t,parts[[3]]$Z_t,tolerance=1e-12)
  expect_equal(parts[[2]]$Z_t[c(3,4)],c(0,0))
})

test_that("six insurance types follow independent future cash-flow ledgers", {
  d <- data.frame(Tx=c(.5,1,1.5,2,2.5,3.5,NA))
  d$Kx <- floor(d$Tx)
  for (type in c("whole","term","deferred","deferred_term","pure_endowment","endowment"))
    for (timing in c("end_of_year","moment_of_death")) for (boundary in c("before_payment","after_payment")) {
      args <- list(.data=d,t=c(0,.5,1,2),i=.04,P=6,premium_unit="annualized",
        benefit=100,type=type,timing=timing,k=4,payment=.25,
        annuity_type="deferred_temporary",h=1,n=2,reserve_timing=boundary,
        not_in_force="zero")
      out <- do.call(mc_reserve,args)
      expected <- lapply(seq_len(nrow(out)),function(j) {
        row <- ceiling(j/4); T <- d$Tx[row]; K <- d$Kx[row]; t <- out$t[j]
        if (is.na(T)) return(c(NA,NA,NA,NA))
        if (T<=t) return(c(0,0,0,0))
        death <- if(timing=="end_of_year") K+1 else T
        event <- switch(type,whole=TRUE,term=death<=2,deferred=death>1,
          deferred_term=death>1&&death<=3,
          pure_endowment=if(timing=="end_of_year") K>=2 else T>2,endowment=TRUE)
        date <- if(type=="pure_endowment") 2 else if(type=="endowment") min(death,2) else death
        eligible <- if(boundary=="before_payment") date>=t else date>t
        Z <- if(event&&eligible) 100*1.04^(-(date-t)) else 0
        # Due payments at 1,1.25,...,2.75, strictly before death.
        dates <- seq(1,2.75,by=.25)
        dates <- dates[dates<T & if(boundary=="before_payment") dates>=t else dates>t]
        Y <- sum(.25*1.04^(-(dates-t)))
        income <- sum(1.5*1.04^(-(dates-t)))
        c(Z,Y,income,Z-income)
      })
      expected <- do.call(rbind,expected)
      expect_equal(out$Z_t,expected[,1],tolerance=1e-12)
      expect_equal(out$Y_t,expected[,2],tolerance=1e-12)
      expect_equal(out$future_pv_premiums,expected[,3],tolerance=1e-12)
      expect_equal(out$L_t,expected[,4],tolerance=1e-12)
    }
})

test_that("premium annuity types retain explicit issue calendars", {
  d <- data.frame(Kx=c(0,1,2),Tx=c(.5,1.5,2.5))
  for(type in c("whole","temporary","deferred","deferred_temporary","certain","guaranteed"))
    for(timing in c("due","immediate")) {
      args <- list(.data=d,t=c(0,.5),i=-.02,P=3,premium_unit="per_payment",
        benefit=0,k=2,payment=7,annuity_type=type,premium_timing=timing,
        n=2,n_guar=1,h=.125,not_in_force="zero")
      out <- do.call(mc_reserve,args)
      expected <- vapply(seq_len(nrow(out)),function(j) {
        T <- d$Tx[ceiling(j/2)]; t <- out$t[j]
        if(T<=t) return(0)
        origin <- if(type %in% c("deferred","deferred_temporary")) .125 else 0
        dates <- origin+(if(timing=="due") 0:12 else 1:12)/2
        alive <- if(type=="certain") rep(TRUE,length(dates)) else
          if(timing=="due") dates<T else dates<=T
        if(type=="guaranteed") alive <- alive | if(timing=="due") dates<1 else dates<=1
        if(type %in% c("temporary","deferred_temporary","certain"))
          alive <- alive & if(timing=="due") dates<origin+2 else dates<=origin+2
        sum(3*.98^(-(dates[alive & dates>=t]-t)))
      },0)
      expect_equal(out$future_pv_premiums,expected,tolerance=1e-12)
    }
})

test_that("control coercions fail before overflow and preserve valid m", {
  d <- data.frame(Kx=1,Tx=1.5)
  expect_error(mc_reserve(d,i=0,P=0,k=3e9),"k|range")
  expect_error(mc_reserve(d,i=0,P=0,k=1e-12),"k|positive")
  expect_error(mc_reserve(d,i=0,P=0,m=1e-12),"m|positive")
  out <- mc_reserve(d,i=.05,i_type="nominal_interest",m=3e9,P=0)
  expect_equal(out$m,3e9)
})

test_that("estimated premiums retain finite monetary income at extreme scales", {
  d <- data.frame(Kx=0,Tx=.5)
  for(amount in c(1e-300,1e300)) {
    basis <- 1/amount
    out <- mc_reserve(d,t=c(0,1),i=0,benefit=amount,payment=basis,not_in_force="zero")
    expect_equal(out$premium_per_payment/amount,c(1,1),tolerance=1e-12)
    expect_equal(out$future_pv_premiums/amount,c(1,0),tolerance=1e-12)
    expect_equal(out$L_t/amount,c(0,0),tolerance=1e-12)
    priced <- mc_premium(data.frame(Kx=0,Tx=.5,pv_benefit=amount,pv_annuity=basis),
                          annuity_payment=basis)
    chained <- mc_reserve(priced,t=0,i=0,benefit=amount,payment=basis)
    expect_equal(chained$future_pv_premiums/amount,1,tolerance=1e-12)
    expect_lt(abs(chained$L_t)/amount,1e-12)
  }
})

test_that("coincident benefit and premium can cancel despite separate overflow", {
  out <- mc_reserve(data.frame(Kx=1,Tx=1),t=c(0,.5),i=-1000,i_type="force",
    P=100,premium_unit="per_payment",benefit=100,timing="moment_of_death",
    premium_timing="immediate")
  expect_true(is.infinite(out$Z_t[1]))
  expect_true(is.infinite(out$future_pv_premiums[1]))
  expect_identical(out$L_t,c(0,0))
})

test_that("in-force masks retain their documented complete and curtate bases", {
  d <- data.frame(Kx=c(1,NA),Tx=c(1.5,NA))
  for(basis in c("complete","curtate")) for(mode in c("na","zero")) {
    out <- mc_reserve(d,t=c(0,1,2),i=0,P=2,benefit=10,in_force_basis=basis,
                      not_in_force=mode)
    expect_identical(out$in_force,c(TRUE,TRUE,FALSE,NA,NA,NA))
    expect_true(all(is.na(out$L_t[4:6])))
    if(mode=="zero") expect_identical(out$L_t[3],0) else expect_true(is.na(out$L_t[3]))
  }
})
