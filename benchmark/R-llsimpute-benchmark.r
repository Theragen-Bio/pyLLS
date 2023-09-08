# This script returns computation time and NRMSE for R-data
# You need to run pyLLS-notebook.ipynb before running this script.

# Library installation
# update R packages
if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install('BiocGenerics')
    BiocManager::install('Rcpp')
}
if (!require("pcaMethods", quietly = TRUE)){
    BiocManager::install("pcaMethods")
}

if (!require("data.table", quietly = TRUE)){
    install.packages('data.table')
}

if (!require("progress", quietly = TRUE)){
    install.packages('rlang')
    install.packages('vctrs')
    install.packages('lifecycle')
    install.packages('pkgconfig')
    install.packages('R6')
    install.packages('crayon')
    install.packages('prettyunits')
    install.packages('progress',dependecies=T)
}

library(pcaMethods)
library(data.table)
library(progress)

# parameter setting
params=c(
    # optimize path
    'tcga-exp'='/data03/project/sjoh/00_tools/LLSimpute/pyLLS/package/pyLLS/benchmark/tcga-exp-benchmark.csv'
    ,'ccle-exp'='/data03/project/sjoh/00_tools/LLSimpute/pyLLS/package/pyLLS/benchmark/ccle-exp-benchmark.csv'
    ,'missing-gene-entry'='/data03/project/sjoh/00_tools/LLSimpute/pyLLS/package/pyLLS/benchmark/missing-gene-list.txt'
    ,'outpath'='/data03/project/sjoh/00_tools/LLSimpute/pyLLS/package/pyLLS/benchmark/'
)

# Load data
indata=list()
indata[['tcga-exp']]=data.table::fread(params['tcga-exp'],data.table=F,sep=',')
indata[['ccle-exp']]=data.table::fread(params['ccle-exp'],data.table=F,sep=',')
indata[['missing-gene-entry']]=data.table::fread(input=params['missing-gene-entry'],sep='\t',header=F)

# preprocessing
parsing_missing_gene=function(x=indata[['missing-gene-entry']]){
    x1=list()
    for(row in 1:nrow(x)){
        x1[[row]]=unlist(strsplit(as.character(x[row,1]),','))
    }
    return(x1)
}
pdata=list()
pdata[['missing-gene-entry']]=parsing_missing_gene(x=indata[['missing-gene-entry']])
# Do computation
# ccle-exp and tcga-exp were already processed by pyLLS-notebook.ipynb.
# x=indata[['ccle-exp']];y=indata[['tcga-exp']];minM=10;maxM=20;by=10;times=5;k=10;metric='correlation'
computation_time_for_rLLS_with_fixedK=function(
        x=indata[['ccle-exp']],y=indata[['tcga-exp']],
        missing_gene=pdata[['missing-gene-entry']],
        k=10,metric='correlation'
    ){
        # This function conducts imputation with ccle-exp and tcga-exp.
        # ccle-exp and tcga-exp will be used as reference and target.
        # minM and maxM = The number of missing genes
        # k = The number of probe genes to be evaluated.
        # prepare missing number range
        n_of_missing_genes = computation_time = c()
        # Prepare data
        rownames(x)=x$V1
        rownames(y)=y$V1
        overlap_gene=intersect(rownames(x),rownames(y))
        x1=x[overlap_gene,2:ncol(x)]
        y1=y[overlap_gene,2:ncol(y)]
        x2=t(cbind(x1,y1)) # pcaMethods::LLSimpute requires sample x gene matrix
        missing_samples=colnames(y1)
        # Save result
        nrmse_list = missing_num_list = time_list =c()
        progress = progress_bar$new(total=length(missing_gene),clear=F,format = " Progress: [:bar] :percent, Estimated completion time: :eta")
        for(miss in missing_gene){
            # Make missing gene table
            x3=x2
            x3[missing_samples,miss]=NA
            # Pearson correlation
            # allVariables = F indicates that Linear regression will be conducted for the complete genes.
            start=Sys.time() # Start time
            estimate = pcaMethods::llsImpute(x3,k=k,correlation='pearson',allVariables = F) # gene in column and sample in rows
            pred = completeObs(estimate)
            end=Sys.time() # End time
            # Calculate NRMSE (normalized root mean squar error)
            answer=round(x2[missing_samples,miss],5)
            pred=round(pred[missing_samples,miss],5)
            nrmse=(mean((answer-pred)^2)^(1/2))/sd(answer)
            # store result
            time_list=c(time_list,as.numeric(end-start))
            nrmse_list=c(nrmse_list,nrmse)
            missing_num_list=c(missing_num_list,length(miss))
            progress$tick()
        }
    df=data.frame(matrix(data=NA,nrow=length(missing_num_list),ncol=3),stringsAsFactors = F)
    colnames(df)=c('missing_num','nrmse','time')
    df['missing_num']=missing_num_list
    df['nrmse']=nrmse_list
    df['time']=time_list       
    return(df)
    }

rLLSresult=computation_time_for_rLLS_with_fixedK(
        x=indata[['ccle-exp']],y=indata[['tcga-exp']],
        missing_gene=pdata[['missing-gene-entry']],
        k=10,metric='correlation'
    )
output_file=paste0(params['outpath'],'/R-llsimpute-test-result.csv')
data.table::fwrite(rLLSresult,file=output_file,sep=',',row.names=F)

#-----------------------------------------------------------------------
# Set k=10, missing gene number = 100, sample-sizes (from 10 to 100 by 50).
#-----------------------------------------------------------------------
computation_time_for_rLLS_with_various_sample_size=function(
        x=indata[['ccle-exp']],y=indata[['tcga-exp']],
        missing_gene=pdata[['missing-gene-entry']],
        k=10,missing_gene_num=100,by=50,min_sample=10,max_sample=100,
        metric='correlation'
    ){
        # This function conducts imputation with ccle-exp and tcga-exp.
        # ccle-exp and tcga-exp will be used as reference and target.
        # min_sample and max_sample = range of the sample sizes sample_sizes=seq(0,max_sample,by)
        # sample_sizes[1]=min_sample
        # k = The number of probe genes to be evaluated.
        # missing_gene_num
        # prepare missing number range
        # Prepare data
        rownames(x)=x$V1
        rownames(y)=y$V1
        overlap_gene=intersect(rownames(x),rownames(y))
        x1=x[overlap_gene,2:ncol(x)]
        y1=y[overlap_gene,2:ncol(y)]
        x2=t(cbind(x1,y1)) # pcaMethods::LLSimpute requires sample x gene matrix
        missing_sample_sizes=seq(0,max_sample,by)
        missing_sample_sizes[1]=min_sample
        missing_samples=colnames(y1)
        # Save result
        missing_sample_list=nrmse_list = time_list =c()
        progress = progress_bar$new(total=length(missing_gene),clear=F,format = " Progress: [:bar] :percent, Estimated completion time: :eta")
        # Select missing_gene_list
        missing_num_list=0
        for(miss in missing_gene){
            if(length(miss)!=missing_gene_num){
                next
            }else{
                # Make missing gene table
            for(miss_size in missing_sample_sizes){
                missing_num_list=missing_num_list+1
                x3=x2
                x3[missing_samples[1:miss_size],miss]=NA
                keep_sample=c(missing_samples[1:miss_size],colnames(x1))
                x3=x3[keep_sample,]
                # Pearson correlation
                # allVariables = F indicates that Linear regression will be conducted for the complete genes.
                start=Sys.time() # Start time
                estimate = pcaMethods::llsImpute(x3,k=k,correlation='pearson',allVariables = F) # gene in column and sample in rows
                pred = completeObs(estimate)
                end=Sys.time() # End time
                # Calculate NRMSE (normalized root mean squar error)
                answer=round(x2[missing_samples[1:miss_size],miss],5)
                pred=round(pred[missing_samples[1:miss_size],miss],5)
                nrmse=(mean((answer-pred)^2)^(1/2))/sd(answer)
                # store result
                time_list=c(time_list,as.numeric(end-start))
                nrmse_list=c(nrmse_list,nrmse)
                missing_sample_list=c(missing_sample_list,miss_size)
            }
            progress$tick()
            }
        }
    df=data.frame(matrix(data=NA,nrow=missing_num_list,ncol=4),stringsAsFactors = F)
    colnames(df)=c('missing_gene_num','nrmse','time','missing_sample_size')
    df['missing_gene_num']=missing_gene_num
    df['nrmse']=nrmse_list
    df['time']=time_list  
    df['missing_sample_size']=missing_sample_list
    return(df)
}

rLLSresult_diff_size=computation_time_for_rLLS_with_various_sample_size(
        k=10,missing_gene_num=10,by=50,min_sample=10,max_sample=100,metric='correlation'
    )
output_file=paste0(params['outpath'],'/R-llsimpute-test-result-diff-sizes.csv')
data.table::fwrite(rLLSresult_diff_size,file=output_file,sep=',',row.names=F)