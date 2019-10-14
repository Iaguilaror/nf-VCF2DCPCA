library(ca)
library(ggplot2)
library("cowplot")
library("dplyr")

data_transformation <- function(D){
  # Assertion: stop if D is not a matrix
  stopifnot(is.matrix(D))
  # Store number of columns, number of rows
  NC = dim(D)[2];NR = dim(D)[1];NRC = NC*NR
  GM = mean(D)
  AR = rowMeans(D);AC = colMeans(D)
  SDR = apply(D,1,sd);SDC = apply(D,2,sd)
  C = D
  C = t(t(D)-AC)-AR+GM

  # ANOVA
  # Compute sum of squares for rows and columns
  DSSR = NC*sum((AR - GM)^2);DSSC = NR*sum((AC - GM)^2);DSST = sum((D-GM)^2)
  ## Sum of squares for RxC interactions
  DSSRXC = DSST-DSSC-DSSR
  ## Degrees of freedom
  DDFT = NRC - 1;DDFR = NR - 1;DDFC = NC - 1;DDFRXC = (NR-1)*(NC-1)
  
  return(list(C = C,DSST=DSST,DSSR=DSSR,DSSC=DSSC,DDFT=DDFT,DDFR=DDFR,DDFC=DDFC, 
              DSSRXC=DSSRXC,DDFRXC=DDFRXC))
}

PCA <- function(C){
  # Need call from previous function instead of re-extract from matrix C.
  NR = dim(C)[1]
  NC = dim(C)[2]
  STEST = 1e-8
  # Initialized matrix PR, PC, and vector EV.
  # Determine the limit for the number of components.
  NXL = min(NR-1, NC-1,7)
  PR = matrix(0, nrow = NR, ncol = NXL)
  PC = matrix(0, nrow =NXL, ncol = NC)
  EV = rep(0, NXL)
  # Store temporary information NX ITER STEST SFT
  temp = matrix(0, ncol = 4, nrow = NXL)
  colnames(temp) = c("NX","ITER","STEST","SFT")
  SKIP = FALSE
  col_sd_flag = FALSE
  for(i in 1:NC){
    if(is.na(sd(C[,i]))){
      col_sd_flag = TRUE
    }
  }
  
  row_sd_flag = FALSE
  for(j in 1:NR){
    if(is.na(sd(C[j,]))){
      row_sd_flag = TRUE
    }
  }
  if(col_sd_flag==TRUE){SKIP=TRUE}
  
  else if(row_sd_flag==TRUE){SKIP=TRUE}
  
  else{
    for(NX in 1:NXL){
      SFT = 1
      ITER = 0
      # Initialize VR with uniform random function
      VR = runif(NR, min = -0.5, max = 0.5)
      SR = VR
      while(ITER<200 & SFT>STEST){
        ITER = ITER + 1
        # Compute VC from VR
        VC = t(VR)%*%C
        EVP = sum(VC^2)
        SVP = sqrt(EVP)
        VC = VC/SVP
        # Recalculate VR from VC by reverse sum
        VR = C%*%t(VC)
        EVP = sum(VR^2)
        SVP = sqrt(EVP)
        VR = VR/SVP
        SFT = max(abs(VR-SR))
        # Retain VR as SR
        SR = VR
      }
      # Store eigenvalue first
      EV[NX] = EVP
      
      SQVP = sqrt(SVP)
      VR = VR*SQVP
      VC = VC*SQVP
      
      PR[,NX] = VR
      PC[NX,] = VC
      
      # Write NX, ITER, STEST, SFT into a matrix
      temp[NX,] = c(NX, ITER, STEST, SFT)
      # Remove NX from data copy
      C = C - VR%*%VC
    }
  }
  return(list(PC = PC,PR = PR, EV = EV, SKIP=SKIP))
}

makeplot<-function(PC,PR,NCOL,NAME_COL,NROW,NAME_ROW,TITLE,NXL){
  #Choose the variant and PCs
  choice = 4
  target.PR = PR[,(abs(NXL*choice)-NXL+1):(NXL*abs(choice))]
  target.PC = PC[,(abs(NXL*choice)-NXL+1):(NXL*abs(choice))]

  return(
    list( all_PCA_values_for_rows=target.PR, all_PCA_values_for_cols=target.PC)
    )
}

write_table <- function(anova_result,pca_result,NCOL,NAME_COL,NROW,NAME_ROW,TITLE,NXL){
  options(scipen=20)
  data = anova_result$C
  result = pca_result

  # Double Centered (AMMI)
    EV=rep(0,NXL)
    df=rep(0,NXL)
    for (i in 1:NXL){
      EV[i]=result$EV[i]
      df[i]=NCOL+NROW-2*i-1
    }
    EV_RES=anova_result$DSSRXC-sum(EV)
    df_res=anova_result$DDFRXC-sum(df)
    INTERACTION=paste(substr(NAME_COL,1,1),"x",substr(NAME_ROW,1,1),sep="")
    Source = c("Total",paste(" ",NAME_COL),paste(" ",NAME_ROW),paste(" ",INTERACTION))
    for(i in 1:NXL){
      chr = paste("  ",paste("IPC",i,sep = ''))
      Source = c(Source,chr)
    }
    Source = c(Source,paste("  ","Residual"))
    
    df<-round(c(anova_result$DDFT,anova_result$DDFC,anova_result$DDFR,anova_result$DDFRXC,df,df_res),5)
    SS<-formatC(as.numeric(c(anova_result$DSST,anova_result$DSSC,anova_result$DSSR,anova_result$DSSRXC,EV,EV_RES)),digits = 5,format = 'f')
    table=format(data.frame(df,SS,row.names = Source),justify="right",digits=11)
    name=paste("ANOVA for",TITLE,"analysis 4: Double Centered (AMMI)")
    cat(name,file= paste(TITLE,"anova6.txt",sep=""),append = F,sep='\n')
    cat("-----------------------------------------------",file= paste(TITLE,"anova6.txt",sep=""),append = TRUE,sep='\n')
    cat(capture.output(table), file = paste(TITLE,"anova6.txt",sep=""),append= TRUE, sep = '\n')
    cat("-----------------------------------------------",file= paste(TITLE,"anova6.txt",sep=""),append = TRUE,sep='\n')

    ## return the anova table for extracting SS values
    ## get the DF with the SS values
    anova_table <- data.frame(df,SS,row.names = Source)
    return(anova_table)
}

plot_request<-function(flag,PC,PR,NCOL,NAME_COL,NROW,NAME_ROW,TITLE,NXL){
  PCA_data_list <- makeplot(PC,PR,NCOL,NAME_COL,NROW,NAME_ROW,TITLE,NXL)
  return(PCA_data_list)
}

PCA7<-function(datamatrix){
  
  ## read the matrix and define some values of interest
  ## originally mymatrix.mat is read as dataframe, to get number of rows and cols
  mymatrix.df <- read.table(file = datamatrix, header = F, sep = "\t")
  ## set useful values from dataframe
  NCOL= ncol(mymatrix.df)
  NAME_COL= "COLN"
  NROW= nrow(mymatrix.df)
  NAME_ROW = "ROWS"
  TITLE= "DATA"
  ## delete dataframe from memory
  rm(mymatrix.df)
  ## scan file as legacy from Gauch
  data_vector <- scan(file = datamatrix)
  # Construct a matrix using the vector
  D = matrix(data_vector, nrow = NROW, ncol = NCOL, byrow = TRUE)
  
  NXL = min(nrow(D)-1, ncol(D)-1,7)
  PR = matrix(NA,nrow=NROW,ncol = 7*NXL)
  PC = matrix(NA,nrow=NCOL,ncol = 7*NXL)
  print("Read data. Analyses begin.")
  # MAIN FUNCTION STARTS HERE
  ## set i = 4 to call for legacy execution of DCPCA
  i = 4
      anova_result = data_transformation(D)
      data = anova_result$C
      
      pca_result = PCA(data)
      pr = pca_result$PR
      pc = pca_result$PC
      skip = pca_result$SKIP
      low = (NXL*i-NXL+1)
      upper = (NXL*i)
      PR[,low:upper]= pr[,1:7]
      PC[,low:upper]= t(pc[1:7,])
      if(skip!=TRUE){
        anova.df <- write_table(anova_result,pca_result,NCOL,NAME_COL,NROW,NAME_ROW,TITLE,NXL)
      }
      if (i==4){
        title= "Double Centered (AMMI)"
      }
      if(skip==FALSE){
        print(paste("Finished analysis",i,"of 7:",title))
      }
      else{
        print(paste("Skipped analysis",i,"of 7:",title))
      }
    
    print("Tables are ready in the working directory.")
    
    # Make Graphs
    flag=TRUE
    PCA_data_list <- plot_request(flag,PC,PR,NCOL,NAME_COL,NROW,NAME_ROW,TITLE,NXL)
    ## append the anova dataframe
    PCA_data_list[["anova"]] <- anova.df
    return(PCA_data_list)
}

## test run
# useful_data <- PCA7(datamatrix="data/mymatrix.txt")
useful_data <- PCA7(datamatrix="data/transposed_matrix.txt")

## save the anova df, and re-read to facilitate data transformation...
write.table(x = useful_data$anova, file = "anovatable.tsv",
            append = F, quote = F,
            sep = "\t", row.names = T, col.names = F)
## re-read anova
anova.df <- read.table(file = "anovatable.tsv", header = F, sep = "\t")
##remove spaces
anova.df[,1] <- gsub(pattern = " ", replacement = "", x = anova.df[,1])
## rename rows
rownames(anova.df) <- anova.df$V1
## rename useful cols
colnames(anova.df) <- c("original_names","df","SS")
## WARNING - validate with cfresno or ehlemus
cumulative_eigenvalues <- sum(anova.df[grepl(pattern = "IPC.", x=rownames(anova.df)),3])
## calculate % of variance -only is valid for IPC columns
anova.df$variance_percentage <- anova.df$SS / cumulative_eigenvalues * 100
  
## plot with ggplot
## separate rows DF
rows_PCA.df <- data.frame(useful_data$all_PCA_values_for_rows)
# # rename cols to IPC
names(rows_PCA.df) <- paste0("IPC",1:ncol(rows_PCA.df))
# tag rows as ROWS data
rows_PCA.df$TAG <- factor("ROWS",levels = c("ROWS","COLS"))

## separate cols DF
cols_PCA.df <- data.frame(useful_data$all_PCA_values_for_cols)
# # rename cols to PC
names(cols_PCA.df) <- paste0("IPC",1:ncol(cols_PCA.df))
# tag rows as COLS data
cols_PCA.df$TAG <- factor("COLS",levels = c("ROWS","COLS"))
 
## bind in a long dataframe with tagged data
plotable.df <- rbind(rows_PCA.df, cols_PCA.df)

## make a function to create and save plots for pairs of IPC
plotIPC <- function(DF, varA, varB, IPCA, IPCB) {
  ## plotting ----
  ## define comparable axis limits
  ## min value across every PCA value calculated in the dataset
  min_limit <- DF %>% select(-TAG) %>% min()
  max_limit <- DF %>% select(-TAG) %>% max()
  
  ## define axis names
  xaxisname <- paste(IPCA,round(varA, digits = 2),"%")
  yaxisname <- paste(IPCB,round(varB, digits = 2),"%")
  
  ## dot size
  dotsize = 1
  
  ## common theme
  mytheme <- theme(axis.text.y = element_text(size = 10),
                   axis.text.x = element_text(size = 10),
                   axis.title.x = element_text(size = 12),
                   axis.title.y = element_text(size = 12),
                   axis.ticks = element_line(size = 0.1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   legend.position = "none")

  imposedplot.p <- ggplot(DF,
                          aes(x = get(IPCA),
                              y = get(IPCB),
                              color = TAG, shape = TAG)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(size = dotsize) +
    scale_x_continuous(name = xaxisname,
                       limits = c(min_limit, max_limit)) +
    scale_y_continuous(name = yaxisname,
                       limits = c(min_limit, max_limit)) +
    mytheme
  
  ## save plot
  # compose plot name
  plot_filename=paste0(IPCA,"_vs_",IPCB,"_imposedplot.png")
  ggsave(
    filename = plot_filename,
    plot = imposedplot.p,
    device = "png",
    width = 7, height = 7 , units = "in",
    dpi = 300)

  ## compose single plots before biplots  
  ## generate rowsplot
  rowsplot.p <- DF %>%
    filter(TAG == "ROWS") %>%
    ggplot(aes(x = get(IPCA),
               y = get(IPCB)) )+
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(color = "#F8766D", size = dotsize, shape = 16) +
    scale_x_continuous(name = xaxisname,
                       limits = c(min_limit, max_limit)) +
    scale_y_continuous(name = yaxisname,
                       limits = c(min_limit, max_limit)) +
    mytheme

  ## generate colsplot
  colsplot.p <- DF %>%
    filter(TAG == "COLS") %>%
    ggplot(aes(x = get(IPCA),
               y = get(IPCB)) )+
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(color = "#00BFC4", size = dotsize, shape = 17) +
    scale_x_continuous(name = xaxisname,
                       limits = c(min_limit, max_limit)) +
    scale_y_continuous(name = yaxisname,
                       limits = c(min_limit, max_limit)) +
    mytheme

  ## generate grid plot
  biplot.p <- plot_grid(rowsplot.p,
                        colsplot.p,
                        nrow = 1, align = "h", axis = "t")

  # compose plot name
  plot_filename=paste0(IPCA,"_vs_",IPCB,"_biplot.png")
  # save plots
  ggsave(filename = plot_filename,
         plot = biplot.p,
         device = "png",
         width = 14, height = 7 , units = "in",
         dpi = 300)
}

## get max number of IPCs
maxIPC <- ncol(plotable.df)-1
## start already plotted index at 0
plotted_index <- 0

## start looping throug posible pairs of IPC
for(pc_number1 in 1:maxIPC ) {
  for(pc_number2 in 1:maxIPC ) {
    ## compose the pair of IPC being evaluated and store it in the index vector
    index <- sort(c(pc_number1,pc_number2))
    index <- paste(index, collapse = "-")
    ## calculate if the pair of IPC in turn has been seen before in the index
    previousplotted_check <- sum(grepl(pattern = index, plotted_index))
    ## if the pair of IPC are the same, OR it has been previously evaluated, skip to the NEXT pair
    if( pc_number1 == pc_number2 | previousplotted_check != 0){
      next
    } else {
      ## print the pair of IPC being plotted
      print(paste("plotting",pc_number1,"vs",pc_number2))
      ## print the variance % by each IPC
      IPCA_var <- anova.df[colnames(plotable.df)[pc_number1],"variance_percentage"]
      IPCB_var <- anova.df[colnames(plotable.df)[pc_number2],"variance_percentage"]
      print(paste(pc_number1,":",IPCA_var))
      print(paste(pc_number2,":",IPCB_var))
      
      ## invoke plotting function
      ## call function to plot IPC
      plotIPC(DF=plotable.df, varA=IPCA_var, varB=IPCB_var,
              IPCA=colnames(plotable.df)[pc_number1],
              IPCB=colnames(plotable.df)[pc_number2])
      
      ## save the plotted pair in the already plotted index
      plotted_index <- c(plotted_index, index)
    }
  }
}

## create a scree plot
scree.df <- anova.df[grepl(pattern = "IPC.", x = anova.df$original_names)
                     ,c("original_names","variance_percentage")]

scree.p <- ggplot(data = scree.df, aes(x = original_names, y = variance_percentage)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,100),
                     breaks = seq(0,100, by = 5),
                     labels = paste(seq(0,100, by = 5),"%") ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_line(size = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(filename = "_screeplot.png",
       plot = scree.p,
       device = "png",
       width = 7, height = 7 , units = "in",
       dpi = 300)
