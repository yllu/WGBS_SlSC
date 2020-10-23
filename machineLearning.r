library(randomForest)
library(party)
library(caret)
library(pROC)
library(sampling)
library(rpart)
library(rpart.plot)

# random forest
df <- read.table("gene.MS.tab", header = T)
df <- df[,-which(names(df) %in% c("geneID"))]
df$label_expr <- gsub("0", "Not Expressed", df$label_expr)
df$label_expr <- gsub("1", "Expressed", df$label_expr)
df$label_expr <- as.factor(df$label_expr)

df.strata <- strata(df,
                    stratanames = c("label_expr"),
                    size = c(rep(min(nrow(df[df$label_expr == "Expressed",]), nrow(df[df$label_expr == "Not Expressed",])), 2)),
                    method = "srswor")
df.select <- df[rownames(df) %in% df.strata$ID_unit,]

index <- sample(2, nrow(df.select), replace = TRUE, prob = c(0.7,0.3))
traindata <- df.select[index == 1,] 
testdata <- df.select[index == 2,]

rf <- randomForest(label_expr ~ .,
                  data = traindata,
                  mtry = 8,
                  ntree = 200,
		  importance = T,
                  proximity = TRUE)

## importance
pdf("varImpPlot.gene.MS.tab.pdf", width = 5, height = 5)
varImpPlot(rf, type = 1)
dev.off()

## ROC
pre <- predict(rf, type = 'prob', newdata = testdata)
modelroc <- roc(testdata$label_expr, as.numeric(pre[,2]))
pdf("roc.gene.MS.tab.pdf", height = 5, width = 5)
plot(modelroc, print.auc = TRUE, auc.polygon = TRUE, grid = c(0.1,0.2),
     grid.col = c("green", "red"), max.auc.polygon = TRUE,
     auc.polygon.col = "skyblue", print.thres = TRUE,
     xlab = "Specificity (1-FPR)",
     ylab = "Sensitivity (TPR)")
dev.off()

# decision tree
df <- read.table("gene.MS.tab", header = T)
df <- df[,-which(names(df) %in% c("geneID"))]
df$label_expr <- gsub("0","Not Expressed",df$label_expr)
df$label_expr <- gsub("1","Expressed",df$label_expr)
df$label_expr <- as.factor(df$label_expr)

df.tree <- rpart(label_expr ~ mgCG + mgCHG + mgCHH,
                 data = df,
                 method = "class",
                 model = T,
                 parm = list(split = "information"),
                 control = rpart.control(minsplit = 20,
                                         minbucket = 6,
                                         cp = 0.01,
                                         maxcompete = 4,
                                         maxsurrogate = 5,
                                         usesurrogate = 2,
                                         xval = 10,
                                         surrogatestyle = 0,
                                         maxdepth = 30))
										 
rpart.plot(df.tree,
           type = 3,
           extra = 2,
           under = F,
           fallen.leaves = F,
           digits = 2,
           varlen = 0,
           faclen = 0,
           roundint = T,
           tweak = 1,
           clip.facs = F,
           clip.right.labs = T,
           snip = F,
           box.palette = "auto",
           shadow.col = 0)
