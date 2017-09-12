\name{StabMap}
\alias{StabMap}
\title{Study of External Prefernce mapping stability}
\usage{
StabMap(Y,X,S,n,axis=c(1,2),formula_lm, formula_gam,
dimredumethod=1,pred.na=FALSE,nbpoints=50)
}
\arguments{
\item{Y}{a data frame in which each row represents a product
and each column represents hedonic scores of a
given consumer for the products.}

\item{X}{a data frame in which each row represents a product and
each column represents evaluations of a given sensory
or physico-chemical attribute for the products.}

\item{S}{a data frame crossing products and panelists evaluations
of a given sensory or physico-chemical attribute for the products.
It's required when choosing Canonical Analysis ad dimension reduction method.}

\item{n}{number of selection of random samples from hedonic data for which maps are compared before and after LOESS smoothing.}

\item{axis}{a length 2 vector specifying the components of dimension reduction method to plot.}

\item{formula_lm}{a formula expression. formula_lm for regression model fitting such as formula_lm = "~I(F1*F1)+I(F2*F2)+F1*F2". This formula can be used for polynomial models, GLM and bayesian model.}

\item{formula_gam}{formula for GAM model fitting such as
formula_gam="~s(F1,k=2)+s(F2,k=2)"}

\item{dimredumethod}{1 is assigned to Pricipal Component Analysis(PCA). If 2, Multiple Factor Analysis(MFA). If 3, Canonical Analysis(CA)}

\item{pred.na}{boolean, if TRUE, predictions outside the set of liking scores are removed, if FALSE, we keep all predections in and out the set of liking scores}

\item{nbpoints}{an integer. The number of points in the grid-point of the discretized plan. By default = 50.}
}
\value{
sum squared distances between preferences in maps for each prediction model before and after LOESS smoothing.}

\description{
Helps comparing the stability of 
two submaps performed from two subsamples of hedonic data, before and after LOESS smoothing through the calculation of an average preference distances between maps as an indicator of stability.
}

\details{
This function introduces a comparison
approach between maps through resampling approach based
 on the calculation of an average squared distance of preferences(ASDP) as an indicator of  maps rapprochment and stability.}

\author{
Ibtihel Rebhi, Rihab Boubakri, Dhafer Malouche 
}

\examples{
\dontrun{
# Example of smoothed version of EPM 
res_stab= StabMap(Y,X,S,n=2,axis=c(1,2),
formula_lm="~I(F1*F1)+I(F2*F2)+F1*F2",
formula_gam="~s(F1,k=2)+s(F2,k=2)",dimredumethod=1, nbpoints=50)
res_stab=melt(res_stab)
colnames(res_stab)=c("prob", "var", "value")
# Visualization of stability results
gr<-ggplot(res_stab,aes(x=prob,y=value,fill=prob))+geom_boxplot()+xlab("Methods")+
ylab("Average Squared Distance of Preferences")
gr<-gr+theme_bw()+theme(legend.position = "none",
    axis.text.x = element_text(angle = 90))
gr
}
}
