\name{PrefMap}
\alias{PrefMap}
\title{External Prefernce mapping tools}
\usage{
PrefMap(Y,X,S, axis=c(1,2),
formula="~I(F1*F1)+I(F2*F2)+F1*F2", dimredumethod=1,
predmodel=1, nbpoints=50,pred.na =FALSE,
graph.pred =FALSE,graph.map =TRUE, graph.map.3D =FALSE)
}
\arguments{
\item{Y}{a data frame in which each row represents a product
and each column represents hedonic scores of a
given consumer for the products.}

\item{X}{a data frame in which each row represents a product and
each column represents evaluations of a given sensory
or physico-chemical attribute for the products.}

\item{S}{a data frame crossing products and panelists evaluations
of a given sensory or physico-chemical attribute for each product. It's required when choosing Canonical Analysis as dimension reduction method.}

\item{axis}{a length 2 vector specifying the components of dimension reduction method to plot.}

\item{formula}{a formula expression to be introduced. If predmodel=1 , it corresponds to four structures:
formula = "~I(F1*F1)+I(F2*F2)+F1*F2" for quadratic model. formula ="~ F1 + F2" for vector model.
formula="~ F1 + F2 + (F1*F1 + F2*F2)" for circular model.
formula="~I(F1*F1)+I(F2*F2)" for elliptic model.
If predmodel= 2, it corresponds to GAM structure, formula="~s(F1,k=2)+s(F2,k=2)".If predmodel= 3, GLM structure, formula="~I(F1*F1)+I(F2*F2)+F1*F2". If predmodel= 4, Bayes structure, formula="~I(F1*F1)+I(F2*F2)+F1*F2"}

\item{dimredumethod}{1 is assignedd to Principal Component Analysis(PCA). If 2, Multiple Factor Analysis(MFA). If 3, Canonical Analysis(CA)}

\item{predmodel}{the type of prediction model. If predmodel = 1 : Polynomial modelisation.
If predmodel = 2 : Generelized Additive models (GAM).
If predmodel = 3 : Generelized Linear Models (GLM).
If predmodel = 4 : Bayesian models.}

\item{nbpoints}{an integer. The number of points in th grid-point of the discretized plan of the map provided with the two first components. By default 50}

\item{pred.na}{boolean, if TRUE, predictions outside the set of liking scores are removed (generally [0:10]), if FALSE, we keep all predictions in and outside the set of liking scores}

\item{graph.pred}{boolean, if TRUE, plot of surface prediction
of all consumers is posted.}

\item{graph.map}{boolean, if TRUE, extrenal preference mapping in 2 dimensions is posted}

\item{graph.map.3D}{boolean, if TRUE, extrenal preference mapping is 3 dimensions is posted}
}
\value{
graph.pred plot of surface predictions

graph.map plot of External Preference Mapping

graph.map.3D plot of External Preference Mapping in 3D

nb.NA if pred.na = TRUE, it returns a list with the number of predictions outside the liking scores space for each consumer prediction surface.

pos.NA if pred.na=TRUE, it returns a list with positions of predictions outside the liking scores space for each consumer prediction surface.

occur.NA if pred.na=TRUE, it returns a data frame with total number of predictions outside the liking scores space in each point of the prediction surface of all consumers.
}
\description{
function that provides External Preference Mapping technique and visualizes predictions of individual preferences in function of PCA components from quantitative descriptive variables by fitting either Polynomial regressions, GAM, GLM or Bayesian models.
}
\details{
The external preference mapping method is a very common way used
in sensory data analysis fields to explain hedonic judgments towards products in function of analytic attributes of these products in order to understand the consumers' sensory profile and more especially which characteristics influence positively or negatively consumers' liking (Danzart,2004). A Principal Component Analysis is performed (PCA) on the descriptive data. Individual consumer responses are then regressed according to the chosen two principal components of PCA. Cards of response surface of preference or surface of prediction for each consumer are displayed and then superposed to obtain the map of predictions and the map of preferences for all consumers. Areas of rejection and preferences are underlined via showing the percentage of consumers. This function refers to the method introduced by M. Danzart. It was implemented aldo in SensoMineR package.
}
\examples{
\dontrun{
# Example of External mapping performed from PCA
# dimension reduction method and quadratic regression model.
 data(hedo_biscuit)
 data(sens_biscuit)
 rownames(hedo_biscuit)=hedo_biscuit[,1]
 Y=hedo_biscuit[,-1]
 juge_senso_bisc<- summaryBy(. ~ produit + juge,
 data=sens_biscuit, FUN=c(mean),keep.names = TRUE,na.rm=TRUE)
 S=juge_senso_bisc[,-3]
 prod_bisc=summaryBy(. ~ produit, data=sens_biscuit,
 FUN=c(mean),keep.names = TRUE,na.rm=TRUE)
 rownames(prod_bisc)= prod_bisc[,1]
 X= prod_bisc[,-c(1:3)]
# Map from QR model
 map_QR= PrefMap(Y,X,S,axis=c(1,2),
 formula="~I(F1*F1)+I(F2*F2)+F1*F2",
 dimredumethod=1, predmodel=1, nbpoints=50,pred.na =FALSE,
 graph.pred =FALSE, graph.map =TRUE, graph.map.3D =FALSE )
# Total number of NA in each point of the grid-map
 nb.QR=map_QR$occur.NA
}
}
\author{
Ibtihel Rebhi, Rihab Boubakri, Dhafer Malouche 
}
\references{
Sebastien Le, Francois Husson (2008).SensoMineR:
A package for sensory data analysis,
 Journal of Sensory Studies Volume 23, Issue 1, Pages 14-25

Danzart M., Sieffermann J.M., Delarue J. (2004). New developments in preference mapping techniques: finding out a consumer optimal product, its sensory profile and the key sensory attributes. 7th Sensometrics Conference, July 27-30, 2004, Davis, CA
}
\seealso{
SmoothMap, StabMap, SensMapUI
}

