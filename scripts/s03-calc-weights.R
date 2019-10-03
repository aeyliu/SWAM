library(MASS)

#September 13, 2017
#this script is to calculate the weights as per the theory derived for the paper
#the result is as follows:
#for each gene
#must first orthogonalize residual matrix such that E[(Y-Yi)(Y-Yj)] = 0 for all tissues
#weights = E[(Y-Y2)^2...(Y-Yn)^2)] 

args = commandArgs(trailingOnly=TRUE)
#args[1] is the location of the measured expression
#args[2] is the directory of predictions
#args[3] is the target name
#args[4] is the extracted model information
#args[5] is index file
#args[6] is output file location
#args[7] is output for ensg map

#args = vector(length=5)
#args[1] = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/example/Cells_EBV-transformed_lymphocytes_Analysis.expr.txt"
#args[2] = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/SWAM/test/lcl-test/intermediate/"
#args[3] = "TW_Cells_EBV-transformed_lymphocytes_0.5_1KG"
#args[4] = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/SWAM/test/lcl-test/info/"
#args[5] = "/net/snowwhite/home/aeyliu/pima/prediXcan/git/SWAM/test/lcl-test/index.txt"

inv.norm = function(z)
{
  q.z = (rank(z)-0.5)/length(z)
  z.n = qnorm(q.z)
  return(z.n)
}



##############################################################
#process predicted expression

tissue.names = as.character(read.table(args[5],header=F)$V1)

predicted.expression = vector( mode="list",length=length(tissue.names) )
predicted.extra = vector( mode="list",length=length(tissue.names) )
genes.list = NULL
weights.all=NULL

genes.temp = NULL
ensg.temp = NULL

for(i in 1:length(tissue.names))
{
 print( paste("Processing ", i, " of ", length(tissue.names),sep="") )
 path.pred = paste(args[2],tissue.names[i],"/predicted_expression.txt",sep="")

 pred.expr = read.table(path.pred,header=TRUE)
 expr.only = pred.expr[,c(-1,-2)]
 path.weights = paste(args[4],tissue.names[i],".extra.dump",sep="")
 weights = read.table(path.weights,header=FALSE,sep="|") #weight is column 3 (V3)
 genes.temp = c(genes.temp,as.character(weights$V2))
 ensg.temp = c(ensg.temp,as.character(weights$V1))
 predicted.extra[[i]] = weights
 weights$V2 = gsub("-",".",as.character(weights$V2))
 overlap = intersect(names(expr.only),as.character(weights$V1))
 weights = weights[as.character(weights$V1)%in%overlap,]
 expr.only = expr.only[,names(expr.only)%in%overlap]
 weights = weights[order(as.character(weights$V1)),]
 expr.only = expr.only[,order(names(expr.only))]
 names(expr.only) = as.character(weights$V2)
 #remove duplicates
 if(length(which(duplicated(weights$V2))==TRUE)>0)
 {
  expr.only = expr.only[,-which(duplicated(weights$V2)==TRUE)]
  weights = weights[-which(duplicated(weights$V2)==TRUE),]
 }
 pred = data.frame(t(expr.only))
 names(pred) = as.character(pred.expr[,1])

 predicted.expression[[i]] = pred
 genes.list = c(genes.list,as.character(weights$V2))
}


#now to compare the actual expression to the predicted expression
in.expression = read.table(args[1],header=T)

#need to process it so it's the same format as the predicted expression
ensg.map = data.frame(ensg.temp,genes.temp)
ensg.map = ensg.map[!duplicated(ensg.map),]
ensg.dupes = which(duplicated(ensg.map[,1]))
rm.index = NULL
#need to resolve cases where there are multiple mappings to gene - take the most common gene mapping
for(i in 1:length(ensg.dupes))
{
 j = which(ensg.map[,1] == ensg.map[ensg.dupes[i],1])
 n.rep = vector(length=length(j))
 for(k in 1:length(j))
 {
  n.rep[k] = length(which(genes.temp == as.character(ensg.map[j[k],2])))
 }
 rm.index = c(rm.index, j[-which.max(n.rep)])
}
if(length(rm.index)>0)
{
 ensg.map=ensg.map[-rm.index,]
}

process.measured = function(measured.expr)
{
 ids.temp = strsplit(as.character(measured.expr$Id),"[.]")
 ids.temp = unlist(ids.temp)
 ids.temp = ids.temp[seq(1,length(ids.temp),by=2)]

 overlap = intersect(measured.expr[,1],ensg.map$ensg)
 measured = measured.expr[which(measured.expr[,1]%in%overlap),]
 ensg = ensg.map[which(ensg.map$ensg%in%overlap),]
 
 measured = measured[order(as.character(measured[,1])),]
 ensg = ensg[order(as.character(ensg$ensg)),]
 measured[,1] = ensg[,2]
 measured = measured[order(measured[,1]),]
 measured.new = aggregate(measured[,2:dim(measured)[2]],list(measured[,1]),mean)
 measured.final = as.matrix(measured.new[,2:dim(measured.new)[2]])
 rownames(measured.final) = as.character(measured.new[,1])
 colnames(measured.final) = gsub("[.]","-",colnames(measured.final))
 return(measured.final)
}

in.final = process.measured(in.expression)


#############
#set target 
target.final = in.final
target.index = which(tissue.names==args[3])


#########################################
cor.to.pv = function(cors,n=465)
{
 val = abs(cors/sqrt((1-cors^2)/(n-2)))
 return(2*(1-pt(val,n-2)))
}

calc.tissue.values = function(expr.pred,expr.target)
{
 expr.pred = expr.pred[,order(colnames(expr.pred))]
 expr.target = expr.target[,order(colnames(expr.target))] #make sure the colnames match
 
 pred.genes = rownames(expr.pred)
 target.genes = rownames(expr.target)
 overlap = intersect(pred.genes,target.genes)
 expr.pred = expr.pred[which(rownames(expr.pred)%in%overlap),]
 extra.genes = rownames(expr.target)[which(!rownames(expr.target)%in%rownames(expr.pred))]
 extra.matrix = matrix(nrow = length(extra.genes),ncol=dim(expr.pred)[2])
 rownames(extra.matrix) = extra.genes
 colnames(extra.matrix) = colnames(expr.pred)
 #expr.pred = t(apply(expr.pred,1,inv.norm))
 
 expr.pred.df = rbind(expr.pred,extra.matrix)
 expr.pred.df = as.matrix(expr.pred.df[order(rownames(expr.pred.df)),])
 #expr.target = expr.target[order(rownames(expr.target)),]
 #expr.target = as.matrix(expr.target)
 #expr.target.int = t(apply(expr.target,1,inv.norm))
 #tissue.residual.matrix = expr.target.int - expr.pred.df

 return(expr.pred.df)
}



#get value of expression for all tissues
tissue.values = vector(mode="list",length=length(tissue.names))
for(i in 1:length(tissue.names))
{
 print(paste("Processing... ", i, " of ", length(tissue.names),sep=""))
 tissue.values[[i]] = calc.tissue.values(predicted.expression[[i]],target.final)
}


##########
#now get the tissue-specific cv prediction accuracy for all tissues

cv.matrix = matrix(nrow=dim(target.final)[1],ncol = length(tissue.names))

for(i in 1:length(tissue.names))
{
 #print(paste("Processing... ", i, " of ", length(tissue.names),sep="")) #no need, it's fast
 genes.temp = rownames(target.final)
 wts.temp = predicted.extra[[i]][,c(2,3)]
 if(sum(duplicated(wts.temp[,1]))>0)
 {
  wts.temp = wts.temp[-which(duplicated(wts.temp[,1])),]
 }

 overlap.temp = intersect(genes.temp,wts.temp$V2)
 extra.temp = genes.temp[which(!genes.temp%in%overlap.temp)]
 dummy.df = data.frame(extra.temp,rep("NA",length(extra.temp)))
 names(dummy.df) = c("V2","V3")
 wts.temp = wts.temp[which(wts.temp$V2%in%overlap.temp),]
 wts.out = rbind(wts.temp,dummy.df)
 wts.out = wts.out[order(as.character(wts.out[,1])),]
 cv.matrix[,i] = as.numeric(as.character(wts.out[,2]))
}
rownames(cv.matrix) = rownames(target.final)

###########################
#calculate weights with inverse covariance matrix


center.values = function(vec)
{
 return(vec-mean(vec))
}


calc.gene.weights2 = function(gene="",values.list,target,index=NULL,diag.val=0,target.index=NULL) 
{
 if(is.null(index)) #if index isn't specified, determine the index by the gene name
 {
  index = which(rownames(values.list[[1]]) == gene)
 }
 
 gene.values = sapply(values.list,function(x) x[index,])

 X.mat = t(gene.values)
 n = dim(X.mat)[2]
 temp = which(apply(X.mat,1,sd)==0)
 tissue.index = which(!is.na(X.mat[,1]))
 X.mat = t(na.omit(X.mat))
 d2 = dim(X.mat)[2]
 
 if(length(temp)>0)
 {
  temp2 = which(tissue.index%in%temp)
  tissue.index=tissue.index[-temp2]
  X.mat = X.mat[,-temp2]
  d2 = d2 - length(temp2)
 }
 Y = target[index,]
 X.mat = as.matrix(X.mat)
 cv.values = cv.matrix[index,tissue.index] #within tissue cross-validated R-squared
 cv.index = NULL
 if(!is.null(target.index))
 {
  cv.index = which(tissue.index==target.index)
 }

 if(d2>0)
 {

  iden = diag(diag.val,d2)
  D = cbind(Y,X.mat)

  #C = apply(D,2,center.values)
  C = apply(D,2,inv.norm)
  C.cov = cor(C)
  C.x = C.cov[-1,-1]+iden
  C.xy = C.cov[1,-1]
  H = ginv(C.x)
  
  #wts = H%*%as.numeric(C.xy)
  cor.xy = as.numeric(cor(Y,X.mat))
  if(!is.null(cv.index))
  {
   cor.xy[cv.index] = sqrt(cv.values[cv.index])
  }
  wts = H%*%cor.xy

  tissue.weights = vector(length=dim(gene.values)[2]) #vector with all the tissue weights, including the NA tissues
  for(i in 1:length(tissue.index))
  {
   temp = wts[i]
   tissue.weights[tissue.index[i]] = temp
  }
 }
 
 if(!d2>0)
 {
  tissue.weights = as.numeric(vector(length=dim(gene.values)[2])) #vector with all the tissue weights, including the NA tissues
 }

 return(tissue.weights)

}


target.weights = matrix(nrow = dim(tissue.values[[1]])[1], ncol = length(tissue.values))

#now do calculate the weights for each gene
for(i in 1:dim(tissue.values[[1]])[1])
{
 if(i %% 1000 == 0)
 {
  print(paste("Processing... ", i, " of ", dim(tissue.values[[1]])[1],sep=""))
 }
 target.weights[i,] = calc.gene.weights2(values.list=tissue.values,target=target.final,index=i,diag.val=3,target.index=target.index)
}

rownames(target.weights) = rownames(tissue.values[[1]])
target.weights = replace(target.weights,target.weights<0,0)

write.table(target.weights,args[6],quote=FALSE,sep="\t")

write.table(ensg.map,args[7],quote=FALSE,sep="\t",row.names=F)


