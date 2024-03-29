library(spatialreg)
data(oldcol, package="spdep")
lw <- spdep::nb2listw(COL.nb)
f <- formula(CRIME ~ INC + HOVAL)
COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, lw, Durbin=TRUE)
p0 <- predict(COL.SLX, newdata=COL.OLD, listw=lw)
COL.SLX.f <- lmSLX(CRIME ~ INC + HOVAL, data=COL.OLD, lw, Durbin=~ INC)
p0f <- predict(COL.SLX.f, newdata=COL.OLD, listw=lw) 
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, Durbin=TRUE)
p1 <- predict(COL.mix.eig)
p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE)
p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
expect_true(isTRUE(all.equal(p2, p3, check.attributes=FALSE)))
COL.mix.eig.f <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, Durbin=~ INC)
p2f <- predict(COL.mix.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS", legacy.mixed = TRUE) 
p3f <- predict(COL.mix.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS", legacy=FALSE, legacy.mixed = TRUE) 
expect_true(isTRUE(all.equal(p2f, p3f, check.attributes=FALSE)))
COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, Durbin=TRUE)
p1 <- predict(COL.SDerr.eig)
p2 <- predict(COL.SDerr.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE)
p3 <- predict(COL.SDerr.eig, newdata=COL.OLD, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE)
expect_true(isTRUE(all.equal(p2, p3, check.attributes=FALSE)))
COL.SDerr.eig.f <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw, Durbin=~ INC)
p2f <- predict(COL.SDerr.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS", legacy.mixed = TRUE) 
p3f <- predict(COL.SDerr.eig.f, newdata=COL.OLD, listw=lw, pred.type = "TS", legacy=FALSE, legacy.mixed = TRUE) 
expect_true(isTRUE(all.equal(p2f, p3f, check.attributes=FALSE)))
set.seed(1)
drops <- match(sample(row.names(COL.OLD), size=9L), row.names(COL.OLD))
COL_40 <- COL.OLD[-drops,]
nd <- COL.OLD[drops,]
lw_40 <- subset(lw, !(attr(lw, "region.id") %in% attr(lw, "region.id")[drops]), zero.policy=TRUE)
nd_COL.SLX <- lmSLX(CRIME ~ INC + HOVAL, data=COL_40, lw_40, Durbin=TRUE, zero.policy=TRUE)
nd_p0 <- predict(nd_COL.SLX, newdata=nd, listw=lw)
nd_COL.SLX.f <- lmSLX(CRIME ~ INC + HOVAL, data=COL_40, lw_40, Durbin=~ INC, zero.policy=TRUE)
nd_p0f <- predict(nd_COL.SLX.f, newdata=nd, listw=lw)
nd_COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL_40, lw_40, Durbin=TRUE, zero.policy=TRUE)
nd_p1 <- predict(nd_COL.mix.eig)
nd_p2 <- predict(nd_COL.mix.eig, newdata=nd, listw=lw, pred.type = "TS", legacy = TRUE, zero.policy=TRUE)
nd_p3 <- predict(nd_COL.mix.eig, newdata=nd, listw=lw, pred.type = "TS", legacy=FALSE, zero.policy=TRUE)
expect_true(isTRUE(all.equal(nd_p2, nd_p3, check.attributes=FALSE)))
nd_COL.mix.eig.f <- lagsarlm(CRIME ~ INC + HOVAL, data=COL_40, lw_40, Durbin=~ INC, zero.policy=TRUE)
expect_error(nd_p2f <- predict(nd_COL.mix.eig.f, newdata=nd, listw=lw, pred.type = "TS", legacy=TRUE)) # https://github.com/r-spatial/spatialreg/issues/38
expect_error(nd_p3f <- predict(nd_COL.mix.eig.f, newdata=nd, listw=lw, pred.type = "TS", legacy=FALSE)) # https://github.com/r-spatial/spatialreg/issues/38
#expect_true(isTRUE(all.equal(nd_p2f, nd_p3f, check.attributes=FALSE)))
nd_COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL_40, lw_40, Durbin=TRUE, zero.policy=TRUE)
nd_p1 <- predict(nd_COL.SDerr.eig)
nd_p2 <- predict(nd_COL.SDerr.eig, newdata=nd, listw=lw, pred.type = "TS",
 legacy.mixed = TRUE, zero.policy=TRUE)
nd_p3 <- predict(nd_COL.SDerr.eig, newdata=nd, listw=lw, pred.type = "TS",
 legacy=FALSE, legacy.mixed = TRUE, zero.policy=TRUE)
expect_true(isTRUE(all.equal(nd_p2, nd_p3, check.attributes=FALSE)))
nd_COL.SDerr.eig.f <- errorsarlm(CRIME ~ INC + HOVAL, data=COL_40, lw_40, Durbin=~ INC, zero.policy=TRUE)
expect_error(nd_p2f <- predict(nd_COL.SDerr.eig.f, newdata=nd, listw=lw, pred.type = "TS", legacy=TRUE, zero.policy=TRUE)) # https://github.com/r-spatial/spatialreg/issues/38
expect_error(nd_p3f <- predict(nd_COL.SDerr.eig.f, newdata=nd, listw=lw, pred.type = "TS", legacy=FALSE, zero.policy=TRUE)) # https://github.com/r-spatial/spatialreg/issues/38
#expect_true(isTRUE(all.equal(nd_p2f, nd_p3f, check.attributes=FALSE)))


