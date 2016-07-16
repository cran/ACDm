.getLLcall <- function(param, dur, exogenousVar = NULL, model, order, mean = mean(dur), distCode = 1,
                       newDay = c(0), returnMu = TRUE, breakPoints = NULL, forceErrExpec = 1,
                       fixedParam = NULL, fixedParamPos = NULL, trace = 0){  
  
  #combines the param and fixedParam into the full param vector if there are fixed parameters:
  if(length(fixedParamPos) != 0)    
    param <- .returnfixedPara(freePara = param, fixedParam = fixedParam, fixedParamPos = fixedParamPos)
  
  distPara <- .seperateStartPara(param, model, distCode, order)$distStartPara
  
  cFunction <- switch(model,
                      ACD = "getLL_ACDcall",
                      LACD1 = "getLL_LACD1call",
                      LACD2 = "getLL_LACD2call",
                      AMACD = "getLL_AMACDcall",
                      ABACD = "getLL_ABACDcall",
                      BACD = "getLL_BACDcall",
                      SNIACD = "getLL_SNIACDcall",
                      LSNIACD = "getLL_logSNIACDcall",
                      stop("model not supported"))
  
  if(length(exogenousVar) == 0){
    
    if(returnMu){
      if(model %in% c("SNIACD", "LSNIACD")){
        temp <- .Call(cFunction,
                      as.double(dur),
                      as.double(param),                     
                      as.integer(order),
                      as.double(mean),
                      as.integer(distCode),
                      as.double(distPara),
                      as.integer(newDay),
                      as.double(breakPoints),
                      as.integer(forceErrExpec), PACKAGE = "ACDm")
      } else {
        temp <- .Call(cFunction,
                      as.double(dur),
                      as.double(param),                     
                      as.integer(order),
                      as.double(mean),
                      as.integer(distCode),
                      as.double(distPara),
                      as.integer(newDay),
                      as.integer(forceErrExpec), PACKAGE = "ACDm")
      }
      
      if(trace != 0) assign("ACDmOptimTrace", c(get("ACDmOptimTrace", envir = ACDmGlobalEnv), param, -temp[[3]]), envir = ACDmGlobalEnv)
      .getLLcall <- list(LL = temp[[3]], mu = temp[[1]], resi = temp[[2]])
      
    } else{
      if(model %in% c("SNIACD", "LSNIACD")){
        LL <- -.Call(cFunction,
                             as.double(dur),
                             as.double(param),                     
                             as.integer(order),
                             as.double(mean),
                             as.integer(distCode),
                             as.double(distPara),
                             as.integer(newDay),
                             as.double(breakPoints),
                             as.integer(forceErrExpec), PACKAGE = "ACDm")[[3]]
      } else {
        LL <- -.Call(cFunction,
                             as.double(dur),
                             as.double(param),                     
                             as.integer(order),
                             as.double(mean),
                             as.integer(distCode),
                             as.double(distPara),
                             as.integer(newDay),
                             as.integer(forceErrExpec), PACKAGE = "ACDm")[[3]]
      }
      
      
      if(trace != 0) assign("ACDmOptimTrace", c(get("ACDmOptimTrace", envir = ACDmGlobalEnv), param, -LL), envir = ACDmGlobalEnv)
      .getLLcall <- LL
      
    }
  } else { #if there are exogenous variables:
    
    cFunction <- paste0(cFunction, "Ex")
    
    if(returnMu){
      if(model %in% c("SNIACD", "LSNIACD")){
        temp <- .Call(cFunction,
                      as.double(dur),
                      as.double(exogenousVar),
                      as.double(param),                     
                      as.integer(order),
                      as.double(mean),
                      as.integer(distCode),
                      as.double(distPara),
                      as.integer(newDay),
                      as.double(breakPoints),
                      as.integer(forceErrExpec), PACKAGE = "ACDm")
      } else {
        temp <- .Call(cFunction,
                      as.double(dur),
                      as.double(exogenousVar),
                      as.double(param),                     
                      as.integer(order),
                      as.double(mean),
                      as.integer(distCode),
                      as.double(distPara),
                      as.integer(newDay),
                      as.integer(forceErrExpec), PACKAGE = "ACDm")
      }
      
      if(trace != 0) assign("ACDmOptimTrace", c(get("ACDmOptimTrace", envir = ACDmGlobalEnv), param, temp[[3]]), envir = ACDmGlobalEnv)
      .getLLcall <- list(LL = -temp[[3]], mu = temp[[1]], resi = temp[[2]])
      
    } else {
      
      if(model %in% c("SNIACD", "LSNIACD")){
        LL <- -.Call(cFunction,
                     as.double(dur),
                     as.double(exogenousVar),
                     as.double(param),                     
                     as.integer(order),
                     as.double(mean),
                     as.integer(distCode),
                     as.double(distPara),
                     as.integer(newDay),
                     as.double(breakPoints),
                     as.integer(forceErrExpec), PACKAGE = "ACDm")[[3]]
      } else {
        LL <- -.Call(cFunction,
                     as.double(dur),
                     as.double(exogenousVar),
                     as.double(param),                     
                     as.integer(order),
                     as.double(mean),
                     as.integer(distCode),
                     as.double(distPara),
                     as.integer(newDay),
                     as.integer(forceErrExpec), PACKAGE = "ACDm")[[3]]
        
      }
      
      if(trace != 0) assign("ACDmOptimTrace", c(get("ACDmOptimTrace", envir = ACDmGlobalEnv), param, -LL), envir = ACDmGlobalEnv)
      .getLLcall <- LL
      
    }
  }
}