#' The survival nomogram R6 class
#' @description
#' create a R6 class for survival nomogram analysis
#' @export
survNomogram <- R6Class(classname = 'survNomogram',
                        list(time = NULL,
                             event = NULL,
                             TimeInterval = NULL,
                             funLabel = NULL,
                             uni_p_thresh = 0.05,
                             fit = NULL,
                             Nomoscore = NULL,
                             cli_name = NULL))
survNomogram$set('public', 'initialize',
                 function(time, event, TimeInterval = c(12, 24, 36), funLabel = c("1-Years Survival Probability",
                                                                     "2-Years Survival Probability",
                                                                     "3-Years Survival Probability"))
                   {
                   self$TimeInterval <- TimeInterval
                   self$funLabel <- funLabel
                   self$time <- time
                   self$event <- event
                 })

survNomogram$set('public', 'run',
                 function(dt, restep = F)
                   {
                   time <- self$time
                   event <- self$event
                   TimeInterval <- self$TimeInterval
                   funLabel <- self$funLabel

                   surv.form0 <- sprintf('Surv(%s, %s)~1', time, event) %>% as.formula()
                   cox.test <- coxphSeries(surv.form0,
                                             data = dt,
                                             vars = colnames(dt)[-c(1:2)])

                   sel.names <- cox.test$Variable[(which(cox.test$Pvalue < self$uni_p_thresh))]


                   dt.1 <- select(dt, c(event, time, sel.names))

                   surv.form <- sprintf('Surv(%s, %s)~.', time, event) %>% as.formula()
                   if(restep)
                   {
                     cox.nom <- coxph(surv.form, data = dt.1) %>% step()
                   }
                   else {
                     cox.nom <- coxph(surv.form, data = dt.1)
                   }

                   coef.nom <- coefficients(cox.nom)
                   Nomoscore <- paste('Nomoscore',
                                      paste(coef.nom, names(coef.nom), sep = '*') %>% paste(collapse = '+'),
                                      sep = '=')
                   cli_name <- names(coef.nom)

                   self$fit <- cox.nom
                   self$Nomoscore <- Nomoscore
                   self$cli_name <- cli_name

                   invisible(self)
                 })

survNomogram$set('public', 'figure',
                 function(dt, outName)
                   {
                   time <- self$time
                   event <- self$event
                   TimeInterval <- self$TimeInterval
                   funLabel <- self$funLabel

                   cox.nom <- self$fit
                   cli_name <- self$cli_name

                   dt.2 <<- dt %>% select(c(event, time, self$cli_name)) %>% as.data.frame()


                   surv.form <- sprintf('Surv(%s, %s)~%s', time, event, paste(self$cli_name,
                                                                              collapse = '+')) %>% as.formula()

                   cox.os.mul <- coxph(surv.form, data = dt.2) %>% publish()

                   cox.os.mul$regressionTable

                   ddist <<- datadist(dt.2); options(datadist='ddist')

                   f <- cph(surv.form, data=dt.2, surv = T)
                   surv <- Survival(f)  # This would also work if f was from cph
                   nom <- nomogram(f, fun=list(function(x) surv(TimeInterval[1], x),
                                               function(x) surv(TimeInterval[2], x),
                                               function(x) surv(TimeInterval[3], x)), funlabel= funLabel,
                                   lp = F, fun.at = seq(0.1, 0.9, by = 0.1))
                   pptx <- read_pptx()
                   pptx <- add_slide(pptx)
                   ph_with(pptx, value = dml(code = {plot(nom)}),
                           location = ph_location_type(type = 'body'))

                   print(pptx, outName)
                 })
survNomogram$set('public', 'predict',
                 function(dt)
                   {
                   cox.nom <- self$fit
                   event <- self$event
                   time <- self$time
                   nomoscore <- predict(cox.nom, newdata = dt)

                   dt.2 <<- dt %>% select(c(event, time, self$cli_name))


                   surv.form <- sprintf('Surv(%s, %s)~%s', time, event, paste(self$cli_name,
                                                                              collapse = '+')) %>% as.formula()

                   cox.mul <- coxph(surv.form, data = dt.2)

                   out <- outSurv_nomo$new()
                   out$Time <- dt[[time]]
                   out$Event <- dt[[event]]
                   out$Score <- nomoscore

                   out$regTable <- publish(cox.mul)$regressionTable

                   out
                 })

#' The output of predict of survNomogram
#' @export
outSurv_nomo <- R6Class(classname = 'outSurv_nomo',
                        inherit = outSurv,
                        public = list(regTable = NULL
                                      ))

