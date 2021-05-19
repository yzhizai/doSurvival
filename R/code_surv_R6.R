#' The R6 object of survival analysis
#' @description
#' The output of the survival analysis
#' @export
outSurv <- R6Class('outSurv', list(
  Time = NULL,
  Event = NULL,
  Score = NULL,
  threshold = NULL,
  cindex = NULL,
  toDf = function()
  {
    tibble(Time = self$Time,
           Event = self$Event,
           Score = self$Score)
  },
#' @description
#' the method to plot km curve
#' @param outName  The output pptx name
#'
#' @return
  kmplot = function(outName)
  {
    dt <- tibble(Time = self$Time,
                 Event = self$Event,
                 Score = self$Score)
    if(is.null(self$threshold))
    {
      cuts <- surv_cutpoint(dt, time = 'Time', event = 'Event',
                            variables = 'Score')

      cuts <- cuts$cutpoint$cutpoint
      self$threshold <- cuts
    }


    dt$Score <- ifelse(dt$Score < self$threshold, 0, 1)
    # KM PFS plot Radscore--------------------------------------------------------------

    survfit.i <- survfit(Surv(Time, Event)~Score, data = dt)

    survplot.i <- ggsurvplot(survfit.i, data = dt, pval = T)

    p.surv <- dml(ggobj = {survplot.i$plot})
    pptx <- read_pptx()
    pptx <- add_slide(pptx)
    ph_with(pptx, value = p.surv,
            location = ph_location_type(typte = 'body'))

    print(pptx, outName)
  },
#' @description The method to calculate the c-index
#' @return NULL
  cidx_calc = function()
  {
    dt <- tibble(Time = self$Time,
                 Event = self$Event,
                 Score = self$Score)
    cox.fit <- coxph(Surv(Time, Event)~Score, data = dt)
    self$cindex <- cox.fit$concordance
  }
))


#' R6 Class Representing a survival radiomics model
#'
#' @description
#' radiomics-based survival model
#'
#' @details
#' Nothing
#' @export
survRadiomics <- R6Class('survRadiomics',
                         list(time = NULL,
                              event = NULL))
survRadiomics$set('public', 'uni_p_thresh', 0.05)
survRadiomics$set('public', 'cor_p_thresh', 0.9)
survRadiomics$set('public', 's_pre', NULL)
survRadiomics$set('public', 'selNames', NULL)
survRadiomics$set('public', 'fit', NULL)
survRadiomics$set('public', 'cv.fit', NULL)
survRadiomics$set('public', 'lambda', NULL)
survRadiomics$set('public', 'initialize',
                  function(time, event, uni_p_thresh = 0.05, cor_p_thresh = 0.9)
                    {
                    self$time <- time
                    self$event <- event
                    self$uni_p_thresh <- uni_p_thresh
                    self$cor_p_thresh <- cor_p_thresh

                    invisible(self)
                  })
survRadiomics$set('public', 'figure',

                  function(outName)
                    {
                    cv.fit <- self$cv.fit
                    fit <- self$fit
                    pptx <- read_pptx()
                    pptx <- add_slide(pptx)

                    p.lasso <- dml(code = {
                      oldpar <- par(mfrow = c(1, 2))
                      plot(cv.fit)
                      plot(fit, s = cv.fit$lambda.min, xvar = 'lambda')
                      abline(v = log(cv.fit$lambda.min), lty = 2)
                      par(oldpar)
                    })

                    ph_with(pptx, value = p.lasso, location = ph_location_type(type = 'body'))

                    pptx <- add_slide(pptx)
                    coefs <- coefficients(fit, s = cv.fit$lambda.min)

                    dt.coef <- tibble(Variable = coefs@Dimnames[[1]][coefs@i + 1],
                                      Coefficients = coefs@x) %>% arrange(desc(Coefficients))

                    Radscore <-paste('Radscore',
                                     paste(dt.coef$Coefficients, dt.coef$Variable, sep = '*') %>% paste(collapse = '+'),
                                     sep = '=')

                    ph_with(pptx, value = Radscore, location = ph_location_type(type = 'body'))
                    # Coefficients plot -------------------------------------------------------

                    pptx <- add_slide(pptx)

                    p.bar <- dml(ggobj = {
                      ggbarplot(data = dt.coef, x = 'Variable', y = 'Coefficients',
                                width = 0.4, fill = 'blue') +
                        coord_flip() + theme_bw()
                    })
                    ph_with(pptx, value = p.bar, location = ph_location_type(type = 'body'))

                    print(pptx, outName)
                  })
survRadiomics$set('public', 'predict',
                  function(dt.radiomics, dt.clinics)
                    {
                    dt.radiomics.1 <- predict(self$s_pre, dt.radiomics)
                    dt.radiomics.2 <- select(dt.radiomics.1, self$selNames)

                    x <- as.matrix(dt.radiomics.2)
                    radscore <- predict(self$fit, newx = x, s = self$lambda) %>% c()

                    res <- outSurv$new()
                    res$Time <- dt.clinics[[self$time]]
                    res$Event <- dt.clinics[[self$event]]
                    res$Score <- radscore

                    res
                  })
survRadiomics$set('public', 'run',
                  function(dt.radiomics, dt.clinics)
                    {
                    s_pre <- preProcess(dt.radiomics,
                                        method = c('medianImpute',
                                                   'center',
                                                   'scale'))

                    dt.radiomics.1 <- predict(s_pre, dt.radiomics)

                    idx_nz <- nearZeroVar(dt.radiomics.1)

                    if(!is_empty(idx_nz))
                    {
                      dt.radiomics.1 <- dt.radiomics.1[, -idx_nz]
                    }

                    dt.temp <- dt.radiomics.1 %>% add_column(Progress = dt.clinics[[self$event]],
                                                             PFS = dt.clinics[[self$time]],
                                                             .before = 1)

                    cox.test <- coxphSeries(Surv(PFS, Progress)~1,
                                            data = dt.temp,
                                            vars = colnames(dt.temp)[-c(1:2)])

                    sel.names <- cox.test$Variable[which(cox.test$Pvalue < self$uni_p_thresh)]

                    rm(dt.temp)

                    dt.radiomics.2 <- select(dt.radiomics.1, sel.names)

                    idx_exc <- findCorrelation(cor(dt.radiomics.2), cutoff = self$cor_p_thresh)
                    idx_in <- setdiff(1:ncol(dt.radiomics.2), idx_exc)

                    dt.radiomics.3 <- dt.radiomics.2[, idx_in]
                    selNames <- colnames(dt.radiomics.3)

                    self$s_pre <- s_pre
                    self$selNames <- selNames

                    x <- as.matrix(dt.radiomics.3)
                    y <- Surv(dt.clinics[[self$time]], dt.clinics[[self$event]])
                    cv.fit <- cv.glmnet(x, y, family = 'cox', nfolds = 10)
                    fit <- glmnet(x, y, family = 'cox')

                    self$fit <- fit
                    self$cv.fit <- cv.fit
                    self$lambda <- cv.fit$lambda.min
                    invisible(self)
                  })
