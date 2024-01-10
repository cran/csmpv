## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4.5) # Set the width and height

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("csmpv")

## ----eval=FALSE---------------------------------------------------------------
#  # Install devtools package if not already installed
#  options(repos = c(CRAN = "https://cloud.r-project.org"))
#  install.packages("devtools")
#  
#  # Install csmpv package from GitHub
#  devtools::install_github("ajiangsfu/csmpv",force = TRUE)
#  # Using force = TRUE will ensure the installation, overriding any existing versions

## ----eval=FALSE---------------------------------------------------------------
#  # Install remotes package if not already installed
#  install.packages("remotes")
#  # Install csmpv package from GitHub
#  remotes::install_github("ajiangsfu/csmpv",force = TRUE)
#  # Using force = TRUE will ensure the installation, overriding any existing versions

## -----------------------------------------------------------------------------
library(csmpv)
data("datlist", package = "csmpv")
tdat = datlist$training
dim(tdat)
vdat = datlist$validation
dim(vdat)

## -----------------------------------------------------------------------------
Xvars = c("highIPI","B.Symptoms","MYC.IHC","BCL2.IHC", "CD10.IHC","BCL6.IHC",
 "MUM1.IHC","Male","AgeOver60", "stage3_4","PS1","LDH.Ratio1",
 "Extranodal1","Bulk10cm","HANS_GCB", "DTI")

## -----------------------------------------------------------------------------
AgeXvars = setdiff(Xvars, "AgeOver60")

## -----------------------------------------------------------------------------
set.seed(12345)

## -----------------------------------------------------------------------------
temp_dir = tempdir()
# setwd(temp_dir) # this only affect this chunk, not for other part
knitr::opts_knit$set(root.dir = temp_dir)

## ----echo=FALSE---------------------------------------------------------------
options(warn = -1) 

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
bconfirm = confirmVars(data = tdat, biomks = Xvars, Y = "DZsig",
                       outfile = "confirmBinary")

## -----------------------------------------------------------------------------
print(bconfirm$fit)
bconfirm$allplot[[2]]

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
cconfirm = confirmVars(data = tdat, biomks = AgeXvars, Y = "Age",
                       outcomeType = "continuous",
                       outfile = "confirmContinuous")

## -----------------------------------------------------------------------------
print(cconfirm$fit)
cconfirm$allplot[[2]]

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
tconfirm = confirmVars(data = tdat, biomks = Xvars,
                       time = "FFP..Years.", event = "Code.FFP",
                       outcomeType = "time-to-event",
                       outfile = "confirmSurvival")

## -----------------------------------------------------------------------------
print(tconfirm$fit)
tconfirm$allplot[[2]]

## -----------------------------------------------------------------------------
bl = LASSO2(data = tdat, biomks = Xvars, Y = "DZsig",
            outfile = "binaryLASSO2")

## -----------------------------------------------------------------------------
bl$coefs

## -----------------------------------------------------------------------------
cl = LASSO2(data = tdat, biomks = AgeXvars,
            outcomeType = "continuous", Y = "Age",
            outfile = "continuousLASSO2")

## -----------------------------------------------------------------------------
cl$coefs

## -----------------------------------------------------------------------------
tl = LASSO2(data = tdat, biomks = Xvars,
            outcomeType = "time-to-event",
            time = "FFP..Years.",event = "Code.FFP",
            outfile = "survivalLASSO2")

## -----------------------------------------------------------------------------
tl$coefs

## -----------------------------------------------------------------------------
b2fit = LASSO2plus(data = tdat, biomks = Xvars, Y = "DZsig",
        outfile = "binaryLASSO2plus")
b2fit$fit$coefficients

## -----------------------------------------------------------------------------
c2fit = LASSO2plus(data = tdat, biomks = AgeXvars,
                   outcomeType = "continuous", Y = "Age",
                   outfile = "continuousLASSO2plus")
c2fit$fit$coefficients

## -----------------------------------------------------------------------------
t2fit = LASSO2plus(data = tdat, biomks = Xvars,
                   outcomeType = "time-to-event",
                   time = "FFP..Years.",event = "Code.FFP",
                   outfile = "survivalLASSO2plus")
t2fit$fit$coefficients

## -----------------------------------------------------------------------------
bfit = LASSO_plus(data = tdat, biomks = Xvars, Y = "DZsig",
                  outfile = "binaryLASSO_plus", topN = 5)
bfit$fit$coefficients

## -----------------------------------------------------------------------------
cfit = LASSO_plus(data = tdat, biomks = AgeXvars,
                  outcomeType = "continuous", Y = "Age",
                  outfile = "continuousLASSO_plus", topN = 5)
cfit$fit$coefficients

## -----------------------------------------------------------------------------
tfit = LASSO_plus(data = tdat, biomks = Xvars,
                  outcomeType = "time-to-event",
                  time = "FFP..Years.",event = "Code.FFP",
                  outfile = "survivalLASSO_plus", topN = 5)
tfit$fit$coefficients

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
blr = LASSO2_reg(data = tdat, biomks = Xvars, Y = "DZsig",
                 outfile = "binaryLASSO2_reg")

## -----------------------------------------------------------------------------
blr$fit$coefficients

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
clr = LASSO2_reg(data = tdat, biomks = AgeXvars,
                 outcomeType = "continuous", Y = "Age",
                 outfile = "continuousLASSO2_reg")

## -----------------------------------------------------------------------------
clr$fit$coefficients

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
tlr = LASSO2_reg(data = tdat, biomks = Xvars,
                 outcomeType = "time-to-event",
                 time = "FFP..Years.",event = "Code.FFP",
                 outfile = "survivalLASSO2_reg")

## -----------------------------------------------------------------------------
tlr$fit$coefficients

## -----------------------------------------------------------------------------
bxfit = XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig",
                    outfile = "binary_XGBoost")
head(bxfit$XGBoost_score)

## -----------------------------------------------------------------------------
cxfit = XGBtraining(data = tdat, biomks = AgeXvars,
                    outcomeType = "continuous", Y = "Age",
                    outfile = "continuous_XGBoost")
head(cxfit$XGBoost_score)

## -----------------------------------------------------------------------------
txfit = XGBtraining(data = tdat, biomks = Xvars,
                    outcomeType = "time-to-event",
                    time = "FFP..Years.",event = "Code.FFP",
                    outfile = "survival_XGBoost")
head(txfit$XGBoost_score)

## -----------------------------------------------------------------------------
blxfit = LASSO2_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig",
                            outfile = "binary_LASSO2_XGBoost")
head(blxfit$XGBoost_score)

## -----------------------------------------------------------------------------
clxfit = LASSO2_XGBtraining(data = tdat, biomks = AgeXvars,
                            outcomeType = "continuous", Y = "Age",
                            outfile = "continuous_LASSO2_XGBoost")
head(clxfit$XGBoost_score)

## -----------------------------------------------------------------------------
tlxfit = LASSO2_XGBtraining(data = tdat, biomks = Xvars,
                            outcomeType = "time-to-event",
                            time = "FFP..Years.",event = "Code.FFP",
                            outfile = "survival_LASSO2_XGBoost")
head(tlxfit$XGBoost_score)

## ----warning=FALSE------------------------------------------------------------
blpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig",
                                 topN = 5,outfile = "binary_LASSO_plus_XGBoost")
head(blpxfit$XGBoost_score)

## ----warning=FALSE------------------------------------------------------------
clpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = AgeXvars,
                                 outcomeType = "continuous", Y = "Age",
                                 topN = 5,outfile = "continuous_LASSO_plus_XGBoost")
head(clpxfit$XGBoost_score)

## ----warning=FALSE------------------------------------------------------------
tlpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = Xvars,
                                 outcomeType = "time-to-event",
                                 time = "FFP..Years.",event = "Code.FFP",
                                 topN = 5,outfile = "survival_LASSO_plus_XGBoost")
head(tlpxfit$XGBoost_score)

## ----warning=FALSE------------------------------------------------------------
bl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig",
                                 outfile = "binary_LASSO2plus_XGBoost")
head(bl2xfit$XGBoost_score)

## ----warning=FALSE------------------------------------------------------------
cl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = AgeXvars,
                                 outcomeType = "continuous", Y = "Age",
                                 outfile = "continuous_LASSO2plus_XGBoost")
head(cl2xfit$XGBoost_score)

## ----warning=FALSE------------------------------------------------------------
tl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = Xvars,
                                 outcomeType = "time-to-event",
                                 time = "FFP..Years.", event = "Code.FFP",
                                 outfile = "survival_LASSO2plus_XGBoost")
head(tl2xfit$XGBoost_score)

## -----------------------------------------------------------------------------
pbl = LASSO2_predict(bl, newdata = vdat, outfile = "pred_LASSO2_binary")
head(pbl)

## -----------------------------------------------------------------------------
pcl = LASSO2_predict(cl, newdata = vdat, outfile = "pred_LASSO2_cont")
head(pbl)

## -----------------------------------------------------------------------------
ptl = LASSO2_predict(tl, newdata = vdat,
                     outfile = "pred_LASSO2_time_to_event")
head(pbl)

## -----------------------------------------------------------------------------
pblr = rms_model(blr$fit, newdata = vdat, outfile = "pred_LASSO2reg_binary")
head(pblr)

## -----------------------------------------------------------------------------
pclr = rms_model(clr$fit, newdata = vdat,
                 outfile = "pred_LASSO2reg_continuous")
head(pclr)

## -----------------------------------------------------------------------------
ptlr = rms_model(tlr$fit, data = tdat, newdata = vdat,
                outfile = "pred_LASSO2reg_time_to_event")
head(ptlr)

## -----------------------------------------------------------------------------
pbfit = rms_model(bfit$fit, newdata = vdat,
                  outfile = "pred_LASSOplus_binary")

## -----------------------------------------------------------------------------
pcfit = rms_model(cfit$fit, newdata = vdat,
                  outfile = "pred_LASSOplus_continuous")

## -----------------------------------------------------------------------------
ptfit = rms_model(tfit$fit, data = tdat, newdata = vdat,
                  outfile = "pred_LASSOplus_time_to_event")

## -----------------------------------------------------------------------------
p2bfit = rms_model(b2fit$fit, newdata = vdat,
                   outfile = "pred_LASSO2plus_binary")

## -----------------------------------------------------------------------------
p2cfit = rms_model(c2fit$fit, newdata = vdat,
                   outfile = "pred_LASSO2plus_continuous")

## -----------------------------------------------------------------------------
p2tfit = rms_model(t2fit$fit, data = tdat, newdata = vdat,
                   outfile = "pred_LASSO2plus_time_to_event")

## -----------------------------------------------------------------------------
pbxfit = XGBtraining_predict(bxfit, newdata = vdat,
                             outfile = "pred_XGBoost_binary")

## -----------------------------------------------------------------------------
pcxfit = XGBtraining_predict(cxfit, newdata = vdat,
                             outfile = "pred_XGBoost_cont")

## -----------------------------------------------------------------------------
ptxfit = XGBtraining_predict(txfit, newdata = vdat,
                             outfile = "pred_XGBoost_time_to_event")

## -----------------------------------------------------------------------------
pblxfit = XGBtraining_predict(blxfit, newdata = vdat,
                              outfile = "pred_LXGBoost_binary")

## -----------------------------------------------------------------------------
pclxfit = XGBtraining_predict(clxfit, newdata = vdat,
                              outfile = "pred_LXGBoost_cont")

## -----------------------------------------------------------------------------
ptlxfit = XGBtraining_predict(tlxfit, newdata = vdat,
                              outfile = "pred_LXGBoost_time_to_event")

## -----------------------------------------------------------------------------
pblpxfit = XGBtraining_predict(blpxfit, newdata = vdat,
                               outfile = "pred_LpXGBoost_binary")

## -----------------------------------------------------------------------------
pclpxfit = XGBtraining_predict(clpxfit, newdata = vdat,
                               outfile = "pred_LpXGBoost_cont")

## -----------------------------------------------------------------------------
ptlpxfit = XGBtraining_predict(tlpxfit, newdata = vdat,
                               outfile = "pred_LpXGBoost_time_to_event")

## -----------------------------------------------------------------------------
pbl2xfit = XGBtraining_predict(bl2xfit, newdata = vdat,
                               outfile = "pred_L2XGBoost_binary")

## -----------------------------------------------------------------------------
pcl2xfit = XGBtraining_predict(cl2xfit, newdata = vdat,
                               outfile = "pred_L2XGBoost_cont")

## -----------------------------------------------------------------------------
ptl2xfit = XGBtraining_predict(tl2xfit, newdata = vdat,
                               outfile = "pred_L2XGBoost_time_to_event")

## -----------------------------------------------------------------------------
vbl = LASSO2_predict(bl, newdata = vdat, newY = TRUE,
                     outfile = "valid_LASSO2_binary")

## -----------------------------------------------------------------------------
vcl = LASSO2_predict(cl, newdata = vdat, newY = TRUE,
                     outfile = "valid_LASSO2_cont")

## -----------------------------------------------------------------------------
vtl = LASSO2_predict(tl, newdata = vdat, newY = TRUE,
               outfile = "valid_LASSO2_time_to_event")

## -----------------------------------------------------------------------------
vblr = rms_model(blr$fit, newdata = vdat, newY = TRUE,
                 outfile = "valid_LASSO2reg_binary")

## -----------------------------------------------------------------------------
vclr = rms_model(clr$fit, newdata = vdat, newY = TRUE,
                 outfile = "valid_LASSO2reg_continuous")

## -----------------------------------------------------------------------------
vtlr = rms_model(tlr$fit, data = tdat, newdata = vdat, newY = TRUE,
                 outfile = "valid_LASSO2reg_time_to_event")

## -----------------------------------------------------------------------------
vbfit = rms_model(bfit$fit, newdata = vdat, newY = TRUE,
                  outfile = "valid_LASSOplus_binary")

## -----------------------------------------------------------------------------
vcfit = rms_model(cfit$fit, newdata = vdat, newY = TRUE,
                  outfile = "valid_LASSOplus_continuous")

## -----------------------------------------------------------------------------
vtfit = rms_model(tfit$fit, data = tdat, newdata = vdat, newY = TRUE,
                  outfile = "valid_LASSOplus_time_to_event")

## -----------------------------------------------------------------------------
v2bfit = rms_model(b2fit$fit, newdata = vdat, newY = TRUE,
                   outfile = "valid_LASSO2plus_binary")

## -----------------------------------------------------------------------------
v2cfit = rms_model(c2fit$fit, newdata = vdat, newY = TRUE,
                   outfile = "valid_LASSO2plus_continuous")

## -----------------------------------------------------------------------------
v2tfit = rms_model(t2fit$fit, data = tdat, newdata = vdat, newY = TRUE,
                   outfile = "valid_LASSO2plus_time_to_event")

## -----------------------------------------------------------------------------
vbxfit = XGBtraining_predict(bxfit, newdata = vdat, newY = TRUE,
                             outfile = "valid_XGBoost_binary")

## -----------------------------------------------------------------------------
vcxfit = XGBtraining_predict(cxfit, newdata = vdat, newY = TRUE,
                             outfile = "valid_XGBoost_cont")

## -----------------------------------------------------------------------------
vtxfit = XGBtraining_predict(txfit, newdata = vdat, newY = TRUE,
                             outfile = "valid_XGBoost_time_to_event")

## -----------------------------------------------------------------------------
vblxfit = XGBtraining_predict(blxfit, newdata = vdat, newY = TRUE,
                              outfile = "valid_LXGBoost_binary")

## -----------------------------------------------------------------------------
vclxfit = XGBtraining_predict(clxfit, newdata = vdat, newY = TRUE,
                              outfile = "valid_LXGBoost_cont")

## -----------------------------------------------------------------------------
vtlxfit = XGBtraining_predict(tlxfit, newdata = vdat, newY = TRUE,
                              outfile = "valid_LXGBoost_time_to_event")

## -----------------------------------------------------------------------------
vblpxfit = XGBtraining_predict(blpxfit, newdata = vdat, newY = TRUE,
                               outfile = "valid_LpXGBoost_binary")

## -----------------------------------------------------------------------------
vclpxfit = XGBtraining_predict(clpxfit, newdata = vdat, newY = TRUE,
                               outfile = "valid_LpXGBoost_cont")

## -----------------------------------------------------------------------------
vtlpxfit = XGBtraining_predict(tlpxfit, newdata = vdat, newY = TRUE,
                               outfile = "valid_LpXGBoost_time_to_event")

## -----------------------------------------------------------------------------
vbl2xfit = XGBtraining_predict(bl2xfit, newdata = vdat, newY = TRUE,
                               outfile = "valid_L2XGBoost_binary")

## -----------------------------------------------------------------------------
vcl2xfit = XGBtraining_predict(cl2xfit, newdata = vdat, newY = TRUE,
                               outfile = "valid_L2XGBoost_cont")

## -----------------------------------------------------------------------------
vtl2xfit = XGBtraining_predict(tl2xfit, newdata = vdat, newY = TRUE,
                               outfile = "valid_L2XGBoost_time_to_event")

## ----eval = FALSE-------------------------------------------------------------
#  modelout = csmpvModelling(tdat = tdat, vdat = vdat,
#                            Ybinary = "DZsig", varsBinary = Xvars,
#                            Ycont = "Age", varsCont = AgeXvars,
#                            time = "FFP..Years.", event = "Code.FFP",
#                            varsSurvival = Xvars,
#                            outfileName= "all_in_one")

## ----warning=FALSE------------------------------------------------------------
DZlassoreg = csmpvModelling(tdat = tdat, vdat = vdat,
                            Ybinary = "DZsig", varsBinary = Xvars,
                            methods = "LASSO2_reg",
                            outfileName= "just_one")

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
xgobj = XGpred(data = tdat, varsIn = Xvars, 
               selection = TRUE,
               time = "FFP..Years.",
               event = "Code.FFP", outfile = "XGpred")

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
tdat$XGpred_class = xgobj$XGpred_prob_class
training_risk_confirm = confirmVars(data = tdat, biomks = "XGpred_class",
                                    time = "FFP..Years.", event = "Code.FFP",
                                    outfile = "training_riskSurvival",
                                    outcomeType = "time-to-event")
training_risk_confirm[[3]]

## -----------------------------------------------------------------------------
xgNew = XGpred_predict(newdat = vdat, XGpredObj = xgobj)

## ----results = 'hide',message=FALSE, warnings=FALSE---------------------------
vdat$XGpred_class = xgNew$XGpred_prob_class
risk_confirm = confirmVars(data = vdat, biomks = "XGpred_class",
                           time = "FFP..Years.", event = "Code.FFP",
                           outfile = "riskSurvival",
                           outcomeType = "time-to-event")
risk_confirm[[3]]

## -----------------------------------------------------------------------------
devtools::session_info()

