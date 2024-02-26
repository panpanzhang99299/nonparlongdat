testthat::test_that("test for time invariant covariates and time invariant auxiliary variables",
                    {
                      test.df <- read.csv(testthat::test_path("data", "mock_data.csv"))
                      
                      ## Available case analysis (for deriving initial values)
                      test.aca.df <-
                        test.df %>% dyplyr::filter(complete.cases(.)) %>% as.data.frame()
                      lme.aca <-
                        nlme::lme(resp ~ covX + covZ, data = test.aca.df, random = ~ 1 |
                                    id)
                      est.aca <-
                        c(as.numeric(fixed.effects(lme.aca)), as.numeric(nlme::VarCorr(lme.aca)[, 2]))
                      
                      ## Data preparation
                      n.sub <- length(unique(test.df$id))
                      m.visit <- sum(test.df$id == 1)
                      data.y <- test.df$resp
                      data.x <- test.df$covX
                      data.z <- test.df$covZ
                      data.aux <- test.df$aux
                      para.ini <- est.aca
                      
                      ## Method implementation
                      est.prop <- X_timeinv-S_timeinv(
                        n.sub,
                        m.visit,
                        data.y,
                        data.x,
                        data.z,
                        data.aux,
                        para.ini,
                        cov.cont = TRUE,
                        aux.cont = TRUE
                      )
                    })    