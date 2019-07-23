## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("cvGEE")
library("geepack")
library("lattice")
library("splines")

## ---- models_Ex1---------------------------------------------------------
pbc2$serCholD <- as.numeric(pbc2$serChol > 210)
gm1 <- geeglm(serCholD ~ ns(year, knots = c(3, 6), Boundary.knots = c(0, 10)) * drug, 
              family = binomial(), data = pbc2, id = id, 
              corstr = "exchangeable")

gm2 <- geeglm(serCholD ~ ns(year, knots = c(2, 5, 7), Boundary.knots = c(0, 10)) * drug, 
              family = binomial(), data = pbc2, id = id, 
              corstr = "exchangeable")

## ---- cv_gee_Ex1---------------------------------------------------------
plot_data <- cv_gee(gm1, return_data = TRUE)
plot_data$non_linear_I <- plot_data$.score
plot_data$non_linear_II <- unlist(cv_gee(gm2))

## ---- plot_Ex1, eval = TRUE, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
xyplot(non_linear_I + non_linear_II ~ year | .rule, data = plot_data, 
       type = "smooth", auto.key = TRUE, layout = c(3, 1),
       scales = list(y = list(relation = "free")),
       xlab = "Follow-up time (years)", ylab = "Scoring Rules")

## ---- models_Ex2---------------------------------------------------------
aids$CD4count <- aids$CD4 * aids$CD4
aids$obstimef <- factor(aids$obstime)

fm1 <- geeglm(CD4count ~ obstimef, family = poisson(), data = aids, 
              id = patient, corstr = "independence")

fm2 <- update(fm1, corstr = "exchangeable")

fm3 <- update(fm1, corstr = "ar1")

## ---- cv_gee_Ex2---------------------------------------------------------
plot_data <- cv_gee(fm1, return_data = TRUE, max_count = 1000)
plot_data$independence <- plot_data$.score
plot_data$exchangeable <- unlist(cv_gee(fm2, max_count = 1000))
plot_data$ar1 <- unlist(cv_gee(fm3, max_count = 1000))

## ---- plot_Ex2, eval = TRUE, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
xyplot(independence + exchangeable + ar1 ~ obstime | .rule, 
       data = plot_data, type = "smooth", auto.key = TRUE, layout = c(3, 1),
       scales = list(y = list(relation = "free")),
       xlab = "Follow-up time (months)", ylab = "Scoring Rules")

