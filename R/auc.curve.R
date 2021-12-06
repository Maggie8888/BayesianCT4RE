
tvauc= function(probs, true_Y){
        # Check if every prediction is identical and if so, return 0.5
        if(length(unique(probs)) == 1L) return(0.5)
        # Convert true_Y to numeric if it's an ordered factor
        if(is(true_Y, "factor")){
                if(is.ordered(true_Y) & length(levels(true_Y)) == 2) actuals <- as.numeric(true_Y) - 1 else stop("true_Y is type factor, but is unordered. Make it an ordered factor.")
        }

        dt <- as.data.frame(cbind(Pred=probs, Actual=true_Y*1L))
        #setorder(dt, -Pred)
        dt <-dt[order(-dt$Pred),]
        CountTrue=aggregate(dt$Actual, by=list(Category=dt$Pred), FUN=sum)
        CountFalse=dim(dt)[1]-aggregate(dt$Actual, by=list(Category=dt$Pred), FUN=sum)
        # Calculate the CumulativeFalsePositiveRate and CumulativeTruePositiveRate
        CumulativeFPR = cumsum(CountFalse)/sum(CountFalse)
        CumulativeTPR = cumsum(CountTrue)/sum(CountTrue)

        # Calculate AUC ROC
        AdditionalArea = c(head(CumulativeFPR, 1) * head(CumulativeTPR, 1)/2,
                                 (tail(CumulativeFPR, -1) - head(CumulativeFPR, -1)) * (head(CumulativeTPR, -1) + (tail(CumulativeTPR, -1) - head(CumulativeTPR, -1))/2))
        auc=1-tail(cumsum(unlist(AdditionalArea)), 1)
        return(auc=auc)
}


