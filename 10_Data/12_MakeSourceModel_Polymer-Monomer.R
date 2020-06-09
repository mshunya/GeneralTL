#------------------------------------------------
#
#Description: 
#
#Author: Shunya Minami
#
#Update records:
#
#------------------------------------------------

##Library
library.list <- c("glmnet", "mxnet", "ranger", "fields", "rBayesianOptimization")

for(i in 1:length(library.list)){
  eval(parse(text=paste("suppressMessages(library(", library.list[i], "))", sep="")))
}
set.seed(1)
mx.set.seed(1)
##Source File
source("./99_Tool.R")

##Function

##Setting
prop.mat <- t(matrix(c(1, 2, "Mono-DielectricConstant",
                       1, 3, "Mono-RefractiveIndex",
                       1, 8, "Mono-HOMOLUMOGap",
                       2, 3, "Poly-BandGap",
                       2, 4, "Poly-DielectricConstant",
                       2, 5, "Poly-RefractiveIndex"
),
3,6))

source.x.area           <- c(2:291)
source.model.name.list  <- c("LM","RF","NN")  #"LM", "RF", "NN"
source.itr.list         <- c(1)
nbo                     <- 100
source.train.ratio      <- 0.8
bc.flg                  <- 1
norm.flg                <- 1

#Parameter List
alpha.list  <- c(0, 1)
lambda.list <- c(0, 1)
ntree.list  <- c(100L, 1000L)
mtry.list   <- c(10L, 100L)

##Main
for (i.prop in c(1:6)){
  data.flg      <- as.numeric(prop.mat[i.prop, 1])
  source.y.area <- as.numeric(prop.mat[i.prop, 2])
  source.y.name <- prop.mat[i.prop, 3]
  if (data.flg == 1){
    source.data <- read.csv("../10_Data/Monomer-all.csv", header=TRUE)
    source.x.area <- c(4357:5237, 10098:11121)
    source.x <- source.data[, source.x.area]
  } else {
    source.data <- read.csv("../10_Data/Polymer-all.csv", header=TRUE)
    source.x.area <- c(4366:5246, 10107:11130)
    source.x <- source.data[, source.x.area]
  }
  source.y <- source.data[, source.y.area]
  
  tmp.id <- !is.na(source.y)
  source.y <- source.y[tmp.id]
  source.x <- source.x[tmp.id,]
  
  for (source.model.name in source.model.name.list){
    for (itr in source.itr.list){
      message(paste(source.y.name, "  ", source.model.name, "   Itr: ", itr, sep=""))
      t <- proc.time()
      
      ##Data
      num.source <- dim(source.x)[1]
      dim.x      <- dim(source.x)[2]
      
      source.train.id <- sample(1:num.source, floor(num.source*source.train.ratio))
      source.x.train  <- as.matrix(source.x[ source.train.id,])
      source.y.train  <- as.vector(source.y[ source.train.id])
      source.x.test   <- as.matrix(source.x[-source.train.id,])
      source.y.test   <- as.vector(source.y[-source.train.id])
      
      bc.list <- list()
      if (bc.flg==1){
        for (i in 1:dim.x){
          bc <- try(powerTransform(source.x.train[,i], family="yjPower"), silent=TRUE)
          if (class(bc) != "try-error"){
            source.x.train[,i] <- yjPower(source.x.train[,i], lambda=bc$lambda)
            source.x.test[,i]  <- yjPower(source.x.test[,i],  lambda=bc$lambda)
          }
          l <- min(source.x.train[,i])
          u <- max(source.x.train[,i])
          if (l==u) {
            source.x.train[,i] <- scale(source.x.train[,i], center=l, scale=1)
            source.x.test[,i]  <- scale(source.x.test[,i],  center=l, scale=1)
          } else {
            source.x.train[,i] <- scale(source.x.train[,i], center=l, scale=(u-l))
            source.x.test[,i]  <- scale(source.x.test[,i],  center=l, scale=(u-l))
          }
          bc.list[[i]] <- list()
          bc.list[[i]]$bc <- bc
          bc.list[[i]]$l  <- l
          bc.list[[i]]$u  <- u
        }
      }
      
      norm.list <- list()
      source.y.train <- scale(source.y.train)
      y.cent <- attributes(source.y.train)$'scaled:center'
      y.scal <- attributes(source.y.train)$'scaled:scale'
      norm.list$y.cent <- y.cent
      norm.list$y.scal <- y.scal
      
      source.y.test <- (source.y.test-y.cent)/y.scal
      
      #Source Model
      bo.id.s <- sample(1:dim(source.x.train)[1], floor(dim(source.x.train)/2))
      source.x.train.tr <- source.x.train[ bo.id.s,]
      source.x.train.cv <- source.x.train[-bo.id.s,]
      source.y.train.tr <- source.y.train[ bo.id.s]
      source.y.train.cv <- source.y.train[-bo.id.s]
      
      if (source.model.name == "LM"){
        bo.lm.s <- function(alpha, lambda){
          tmp.model <- glmnet(x      = source.x.train.tr,
                              y      = source.y.train.tr,
                              alpha  = alpha,
                              lambda = lambda,
                              family = "gaussian")
          pred <- predict(tmp.model, source.x.train.cv)
          score <- -sqrt(mean((pred - source.y.train.cv)^2))
          return(list(Score=score, Pred=pred))
        }
        tmp <- 0
        while(tmp==0){
          opt.lm.s <- try(BayesianOptimization(bo.lm.s,
                                               bounds      = list(alpha=alpha.list, lambda=lambda.list),
                                               init_points = nbo,
                                               n_iter      = 1,
                                               acq         = "ei",
                                               kappa       = 2.576,
                                               eps         = 0,
                                               verbose     = FALSE),
                          silent = TRUE)
          if (class(opt.lm.s)!="try-error"){
            tmp <- 1
          }
        }
        
        alpha.s  <- opt.lm.s$Best_Par[1]
        lambda.s <- opt.lm.s$Best_Par[2]
        
        source.model <- glmnet(x      = source.x.train,
                               y      = source.y.train,
                               alpha  = alpha.s,
                               lambda = lambda.s,
                               family = "gaussian")
        
        source.y.fits.scal <- predict(source.model, source.x.train)
        source.y.pred.scal <- predict(source.model, source.x.test)
        
        source.param <- paste("alpha=", round(alpha.s, 3), " ,lambda=", round(lambda.s, 3), sep="")
        
        #Importance
        sds <- apply(source.x.train, 2, sd)
        source.importance <- source.model$beta * sds
        
      } else if (source.model.name == "RF"){
        df.train.s    <- data.frame(source.x.train, y=source.y.train)
        df.train.s.tr <- df.train.s[ bo.id.s,]
        df.train.s.cv <- df.train.s[-bo.id.s,]
        
        bo.rf.s   <- function(ntree, mtry){
          tmp.model <- ranger(y~.,
                              data     = df.train.s.tr,
                              num.tree = ntree,
                              mtry     = mtry)
          pred  <- predict(tmp.model, df.train.s.cv)$predict
          score <- -sqrt(mean((pred - source.y.train.cv)^2))
          return(list(Score=score, Pred=pred))
        }
        tmp <- 0
        while(tmp==0){
          opt.rf.s <- try(BayesianOptimization(bo.rf.s,
                                               bounds      = list(ntree=ntree.list, mtry=mtry.list),
                                               init_points = nbo,
                                               n_iter      = 1,
                                               acq         = "ei",
                                               kappa       = 2.576,
                                               eps         = 0,
                                               verbose     = FALSE),
                          silent = TRUE)
          if (class(opt.rf.s)!="try-error"){
            tmp <- 1
          }
        }
        
        ntree.s <- opt.rf.s$Best_Par[1]
        mtry.s  <- opt.rf.s$Best_Par[2]
        
        source.model <- ranger(y~.,
                               data     = df.train.s,
                               num.tree = ntree.s,
                               mtry     = mtry.s,
                               importance = "permutation")
        
        source.param <- paste("ntree=", round(ntree.s, 3), " ,mtry=", round(mtry.s, 3), sep="")
        
        source.y.fits.scal <- predict(source.model, df.train.s)$predict
        source.y.pred.scal <- predict(source.model, as.data.frame(source.x.test))$predict
        
        #Importance
        source.importance <- importance(source.model)
        
      } else if (source.model.name == "NN"){
        #num.layer <- sample(3:5, 1)
        num.layer  <- 3
        num.hidden <- c(rep(1, (num.layer+1)))
        tmp.n      <- dim.x
        for (nl in 1:num.layer){
          tmp.n          <- max(min(max(2,ceiling(runif(1,min=0.8,max=0.9)*tmp.n)),1000),10)
          num.hidden[nl] <- tmp.n
        }
        num.hidden <- c(1000, 200, 50, 1)
        data <- mx.symbol.Variable("data")
        fc   <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=num.hidden[1])
        for (f in 2:(num.layer+1)){
          name.fc  <- paste("fc", f, sep="")
          name.act <- paste("act", (f-1), sep="")
          act      <- mx.symbol.Activation(fc, name=name.act, act_type="relu")
          fc       <- mx.symbol.FullyConnected(act, name=name.fc, num_hidden=num.hidden[f])
        }
        lro <- mx.symbol.LinearRegressionOutput(fc)
        
        source.model <- mx.model.FeedForward.create(lro,
                                                    X                = source.x.train,
                                                    y                = as.vector(source.y.train),
                                                    optimizer        = "sgd",
                                                    ctx              = mx.cpu(),
                                                    num.round        = 150,
                                                    array.batch.size = 16,
                                                    initializer      = mx.init.uniform(0.1),
                                                    learning.rate    = 0.001,
                                                    momentum         = 0.9,
                                                    array.layout     = "rowmajor",
                                                    eval.metric      = mx.metric.rmse,
                                                    verbose          = TRUE)
        
        source.y.fits.scal <- as.vector(predict(source.model, source.x.train, array.layout="rowmajor"))
        source.y.pred.scal <- as.vector(predict(source.model, source.x.test, array.layout="rowmajor"))
        
        source.param <- paste(num.hidden[1],sep="")
        for (p in 2:(num.layer+1)){
          source.param <- paste(source.param, "-", num.hidden[p], sep="")
        }
        
        source.importance <- NA
      }
      
      source.y.fits  <- source.y.fits.scal * y.scal + y.cent
      source.y.pred  <- source.y.pred.scal * y.scal + y.cent
      source.y.train <- source.y.train * y.scal + y.cent
      source.y.test  <- source.y.test  * y.scal + y.cent
      
      #Save
      if (source.model.name == "NN"){
        save.name.source   <- paste("../30_Output/10_Model/00_NN/10_SourceModel_", 
                                    source.y.name, "_", source.model.name, "_itr", itr, sep="")
        mx.model.save(source.model, save.name.source,   0)
      }
      
      model <- list()
      model$train.id   <- source.train.id
      model$bc.flg     <- bc.flg
      model$bc.list    <- bc.list
      model$norm.flg   <- norm.flg
      model$norm.list  <- norm.list
      model$model.name <- source.model.name
      model$model      <- source.model
      model$param      <- source.param
      model$importance <- source.importance
      model$accuracy   <- make.result(source.y.test, source.y.pred)
      model$proc.time  <- (proc.time()-t)[3]
      
      save.name <- paste("../30_Output/10_Model/10_SourceModel_",
                         source.y.name, "_", source.model.name, "_itr", itr, ".RData", sep="")
      save(model, file=save.name)
      
      #Plot
      plot.name <- paste("../30_Output/20_Plot/10_SourceModel_Plot_", 
                         source.y.name, "_", source.model.name, "_itr", itr, ".jpeg", sep="")
      
      jpeg(plot.name, width=800, height=800)
      par(oma=c(0,0,0,0))
      
      xylim.source <- range(source.y.fits, source.y.train, source.y.pred, source.y.test, finite=TRUE)
      grid.plot.s <- seq(-100, 100, by=2)
      grid.plot <- seq(-100, 100, by=0.5)
      
      par(pch=20, cex.main=1.5, font.main=2, cex.axis=1.2, font.axis=2, cex.lab=1.3, font.lab=2)
      
      #Source Model
      plot(x    = source.y.fits,
           y    = source.y.train,
           xlim = xylim.source,
           ylim = xylim.source,
           col  = rgb(1,0,0, alpha=0.2),
           xlab = "",
           ylab = "",
           xaxt="n",
           yaxt="n",
           main = "")
      par(new=TRUE)
      plot(x    = source.y.pred,
           y    = source.y.test,
           xlab = "Prediction",
           ylab = "Observation",
           xlim = xylim.source,
           ylim = xylim.source,
           pch  = 20,
           col  = rgb(0,0,1, alpha=0.2),
           main = "Source Model")
      legend("topleft", legend=make.legend(source.y.test, source.y.pred), cex=1, bg="white")
      abline(a=0, b=1, lty=2)
      #abline(v=grid.plot.s, lty="dotted", col=rgb(0,0,0,alpha=0.2))
      #abline(h=grid.plot.s, lty="dotted", col=rgb(0,0,0,alpha=0.2))
      mtext(text=paste("NN : ", source.param, sep=""),adj=0, outer=FALSE, cex=0.9)
      dev.off()
    }
  }
}