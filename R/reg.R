#' Regression Matrix
#' @description Computes regression such that one variable is regressed over all other variables
#'
#' @param data A dataset
#'
#' @param family Error distribution to be used in the regression model.
#' Defaults to \code{"logistic"}.
#' Set to any family used in function \code{\link[stats]{family}}
#'
#' @param symmetric Should matrix be symmetric?
#' Defaults to \code{TRUE}, taking the mean of the two edge weights
#' (i.e., \code{[i,j]} and \code{[j,i]})
#' Set to \code{FALSE} for asymmetric weights
#' (i.e., \code{[i,j]} does not equal \code{[j,i]})
#'
#' @return A matrix of fully regressed coefficients where
#' one variable is regressed over all others
#'
#' @examples
#' #binarize responses
#' psyb <- ifelse(neoOpen>=4, 1, 0)
#'
#' \dontrun{
#' #perform logistic regression
#' mat <- reg(psyb)}
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @importFrom stats glm
#'
#' @export
#Regression----
reg <- function (data,
                 family = c("binomial" ,"gaussian", "Gamma", "poisson"),
                 symmetric = TRUE)
{
    if(missing(family))
    {family<-"binomial"
    }else(family<-match.arg(family))

    n<-ncol(data)

    data<-as.data.frame(data)

    mat<-matrix(0,nrow=(n-1),ncol=n)

    pb <- txtProgressBar(max=n, style = 3)
    for(i in 1:ncol(data))
    {
        res<-cbind(data[,i],data[,-i])
        mat[,i]<-glm(res,family=family)$coefficients[2:(ncol(data))]

        setTxtProgressBar(pb, i)
    }
    close(pb)

    nmat<-matrix(0,nrow=n,ncol=n)

    for(i in 1:n)
    {
        if(i==1)
        {nmat[,i]<-c(0,mat[,i])
        }else if(i!=n)
        {nmat[,i]<-c(mat[1:(i-1),i],0,mat[i:nrow(mat),i])
        }else if(i==n)
        {nmat[,i]<-c(mat[,i],0)}
    }

    if(symmetric)
    {nmat<-(nmat+t(nmat))/2}

    row.names(nmat)<-colnames(data)
    colnames(nmat)<-colnames(data)

    return(nmat)
}
#----