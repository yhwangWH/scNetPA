# Using R for visualization
library(igraph)
library(ggplot2)
library(ggsci)
library(randomForest)
library(e1071)

edge_list = read.csv('filtered_edge_list.csv')
tf_list = read.csv('tf_list.csv')
node_id_mapping = read.csv('gene_name_id_mapping.csv') # from SCENIC
gene_name = read.table('normalizedDataGeneNameList.txt')
current_gene_exp = read.table('normalizedData.txt')[,1] # 36601


node_id_mapping = node_id_mapping[,c(2,1)] 

node_id_mapping$category = c('gene','TF gene')[1+(node_id_mapping[,2] %in% tf_list[,1])] 
names(node_id_mapping) = c("id","name","category")
edges = edge_list[,2:3]
g = graph.data.frame(edges,directed=T,vertices=node_id_mapping) 
E(g)$weights = edge_list[,4] 


exp_component = g

# exp_component = delete.vertices(exp_component,remove_gene_name) 
Isolated = which(degree(exp_component)==0) 
exp_component = delete.vertices(exp_component, Isolated) 
Isolated = which(degree(exp_component)<2) 
g2 = delete.vertices(exp_component, Isolated)


# extract the largest component
components = decompose.graph(g2) 
vcounts = sapply(components,gorder) 

# take the first component for example
components = components[order(vcounts,decreasing=T)]

example_graph = simplify(components[[1]],,edge.attr.comb='min')
plot(example_graph
     ,edge.arrow.size=0.3 
     ,vertex.size=6.5
     ,vertex.label.cex=0.5
     ,vertex.frame.color=NA
     #,vertex.color=(V(example_graph)$category=='gene')+1
     ,edge.color = c('salmon','grey')[(E(example_graph)$weights>0)+1]
     ,layout=layout_with_lgl(example_graph) )

############################################
############################################
## GO analysis ## 
V(example_graph)$name
GO.barplot=function(file,index,N,main,ylim){
  require(ggplot2)
  
  GO = read.csv(file,header=F)
  
  GO = subset(GO,GO[,2]<0.05)
  
  if( N == 0) {N = dim(GO)[1]}
  GO = GO[1:N,]
  col=rep("black",N)
  col[index] = "black"
  GO = as.data.frame(GO)
  names(GO) = c("GO_name","p_value")
  GO$p_value = -log10(GO$p_value)
  GO$group = col
  fa = factor(GO$GO_name,levels = GO$GO_name)
  GO$GO_name=fa
  rhg_cols = c("#3399ff","#ff6600")
  
  ggplot(GO,aes(x=GO_name,y=p_value,fill=group))+geom_bar(,stat="identity")+scale_fill_manual(values = rhg_cols)+
    theme(axis.title=element_text(size=20,face="bold"),axis.text.y = element_text(size=20),legend.position="none",axis.text.x = element_text(angle = 80, hjust = 1,size=8,color=GO$group))+labs(y="-log(p_value)")+ggtitle(main)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 40, simplify = FALSE), paste, collapse="\n"))+
    ylim(ylim)+geom_hline(aes(yintercept=-log10(0.05)), color="red", linetype="dashed")
}

GO.barplot('GONormed.csv',,30,'GO analysis',c(0,25)) # 
GO.barplot('KEGGNormed.csv',,16,'KEGG analysis',c(0,4)) # 

############################################
############################################
## PPI analysis ##
edge_list = as_edgelist(example_graph)
edge_list  = subset(edge_list,E(example_graph)$weights>0) 
ppi = read.table("9606.ppi.txt") # schema: p1, p2, experiments score, combined score

condense_two_vec = function(vec1,vec2,mode='b'){
  # mode: 
  # b: return 1 if either vec1 or vec2 is not NA
  # n: return the largest element of vec1 and vec2
  res_mat = rbind(vec1,vec2)
  if(mode=='b'){
    res = apply(res_mat,2,function(x){ # 2 indicates cols
      any(!is.na(x))
    })
  }
  if(mode=='n'){
    res = apply(res_mat,2,function(x){
      if(all(is.na(x))){
        return(0)
      }else{
        # na.rm	a logical indicating whether missing values should be removed.
        max(x,na.rm=T)
      }
    })
  }
  return(res)
}

ppi_fraction = function(threshold,ppi,edge_list){
  # threshold: threshold of the experiments score to remove the low quality ppi
  edge_list = as.data.frame(edge_list)
  ppi = as.data.frame(ppi)
  sub_ppi = subset(ppi,ppi[,4]>threshold)
  id_1 = prodlim::row.match(edge_list,sub_ppi[,1:2])
  score_1 = sub_ppi[id_1,4]
  id_2 = prodlim::row.match(edge_list,sub_ppi[,2:1])
  score_2 = sub_ppi[id_2,4]
  binarized_res = condense_two_vec(id_1,id_2,'b')
  score_res = condense_two_vec(score_1,score_2,'n')
  return(list(binarized_res,score_res))
}

background_ppi_fraction = function(threshold,ppi,edge_list){
  # shuffle the paired nodes to generate the background edges
  edge_list[,2] = sample(edge_list[,2])
  return(ppi_fraction(threshold,ppi,edge_list))
}

baseline_ppi_fraction = function(threshold,ppi,edge_list,n=1){
  # randomly select same number of genes from all ppi pool
  # to generate random edges
  n_edges = dim(edge_list)[1]
  sub_ppi = subset(ppi,ppi[,4]>threshold)
  ppi_gene_pool = unique(unlist(sub_ppi[,1:2]))
  res = list()
  for(i in 1:n){
    print(i)
    left_gene = sample(ppi_gene_pool,n_edges) 
    right_gene = sample(ppi_gene_pool,n_edges)
    edge_list = data.frame(left_gene,right_gene)
    res[[i]] = ppi_fraction(threshold,ppi,edge_list)
  }
  return(res)
}


se = function(x){sd(x)/sqrt(length(x))}

fg = ppi_fraction(100,ppi,edge_list) 

bg = list()
for(i in 1:100){
  print(i)
  bg[[i]] = background_ppi_fraction(100,ppi,edge_list)
}


baseline = baseline_ppi_fraction(100,ppi,edge_list,100)

bg_frac_mean = mean(sapply(bg,function(x){
  mean(x[[1]])
}))

bg_frac_se = se(sapply(bg,function(x){
  mean(x[[1]])
}))

baseline_mean = mean(sapply(baseline,function(x){
  mean(x[[1]])
}))

baseline_se = se(sapply(baseline,function(x){
  mean(x[[1]])
}))

df = data.frame(
  gs=c('baseline','random edge','network edge'),    
  mean=c(baseline_mean,bg_frac_mean,mean(fg[[1]])),
  sd=c(baseline_se,bg_frac_se,0)
)


mypal = rev(pal_futurama(alpha=1)(9)[c(2,5,9)]) 

df$'gs' <- factor(as.character(df$'gs'),levels=c('baseline','random edge','network edge'))

p <- ggplot(df, aes(x=gs,y=mean,fill=gs)) + 
  geom_bar(position="dodge", stat="identity")+
  ylab("PPI validated fraction")+
  xlab("Method")+

  geom_errorbar(position=position_dodge(.9),aes(ymin=mean-sd, ymax=mean+sd), width=.3)+ 

  scale_fill_manual(values=mypal)+
 
  theme(text = element_text(size=20))+
  theme_classic()+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(text = element_text(size = 20))                                                                                                           
print(p)
# save.image() 

#####################################
### 
### Predict gene expression from neighbors ###
exp_mat_path = "GSM5687481_k562_rep1_matrix.mtx.gz"
gene_name_path = "GSM5687481_k562_rep1_features.tsv.gz"
exp_mat = Matrix::readMM(exp_mat_path)
exp_mat = t(as.matrix(exp_mat))

gene_name_list = read.table(gene_name_path)

estimate_gene_exp = function(target_gene_name, var_gene_name, gene_name_list, input_exp_mat, estimator){

  target_gene_id = which(gene_name_list[,2] == target_gene_name)
  
  var_gene_id = match(var_gene_name, gene_name_list[,2])
  y = input_exp_mat[,target_gene_id]
  x = input_exp_mat[,var_gene_id]
  df = data.frame(cbind(y,x))
  colnames(df)[-1] = var_gene_name
  fitted_estimator = estimator(y~.,data=df)
  return(fitted_estimator)
}

get_var_gene_name = function(target_gene_name, edge_list){
  return(subset(edge_list[,1],edge_list[,2]==target_gene_name))
}

eval_estimator = function(target_gene_name, var_gene_name, gene_name_list, input_eval_dataset, estimator, log_tranformed=F){

  target_gene_id = which(gene_name_list[,2] == target_gene_name)
  var_gene_id = match(var_gene_name, gene_name_list[,2])

  y = input_eval_dataset[,target_gene_id]
  x = as.data.frame(input_eval_dataset[,var_gene_id])
  colnames(x) = var_gene_name
 
  prediction = predict(estimator, x)
  if(!log_tranformed){
    return( mean(abs(y-prediction)) )
  }else{
    return( mean(abs(exp(y)-exp(prediction))) )
  }
}


in_degree = degree(
  example_graph,
  v = V(example_graph),
  # “in” for in-degree
  mode = c("in"),
  loops = TRUE,
  normalized = FALSE
) 


edge_list = get.edgelist(example_graph) 
target_gene_names = names(which(in_degree>0)) 
var_gene_list = lapply(target_gene_names,function(x){get_var_gene_name(x,edge_list)})
names(var_gene_list) = target_gene_names

# test different algorithms here to predict the gene exp and compare performance
# lm vs svm vs rf
train_dataset = exp_mat[1:2250,] 
eval_dataset = exp_mat[-c(1:2250),]
lm_estimator = list()
rf_estimator = list()
svm_estimator = list()
eval_mae = c()

for(target_gene_name in target_gene_names){ 
  var_gene_name = var_gene_list[[target_gene_name]]
  lm_estimator[[target_gene_name]] = estimate_gene_exp(target_gene_name, var_gene_name, gene_name_list, train_dataset, lm)
  rf_estimator[[target_gene_name]] = estimate_gene_exp(target_gene_name, var_gene_name, gene_name_list, train_dataset, randomForest)
  svm_estimator[[target_gene_name]] = estimate_gene_exp(target_gene_name, var_gene_name, gene_name_list, train_dataset, svm)
  lm_mae = eval_estimator(target_gene_name, var_gene_name, gene_name_list, eval_dataset, lm_estimator[[target_gene_name]])
  rf_mae = eval_estimator(target_gene_name, var_gene_name, gene_name_list, eval_dataset, rf_estimator[[target_gene_name]])
  svm_mae = eval_estimator(target_gene_name, var_gene_name, gene_name_list, eval_dataset, svm_estimator[[target_gene_name]])
  eval_mae = rbind(eval_mae,c(lm_mae,rf_mae,svm_mae))
}

# checkout if the log-transformation is needed
log_lm_estimator = list()
log_rf_estimator = list()
log_svm_estimator = list()
log_eval_mae = c()
log_train_dataset = log(1+train_dataset)
log_eval_dataset = log(1+eval_dataset)
for(target_gene_name in target_gene_names){
  var_gene_name = var_gene_list[[target_gene_name]]
  log_lm_estimator[[target_gene_name]] = estimate_gene_exp(target_gene_name, var_gene_name, gene_name_list, log_train_dataset, lm)
  log_rf_estimator[[target_gene_name]] = estimate_gene_exp(target_gene_name, var_gene_name, gene_name_list, log_train_dataset, randomForest)
  log_svm_estimator[[target_gene_name]] = estimate_gene_exp(target_gene_name, var_gene_name, gene_name_list, log_train_dataset, svm)
  lm_mae = eval_estimator(target_gene_name, var_gene_name, gene_name_list, log_eval_dataset, log_lm_estimator[[target_gene_name]], log_tranformed=T)
  rf_mae = eval_estimator(target_gene_name, var_gene_name, gene_name_list, log_eval_dataset, log_rf_estimator[[target_gene_name]], log_tranformed=T)
  svm_mae = eval_estimator(target_gene_name, var_gene_name, gene_name_list, log_eval_dataset, log_svm_estimator[[target_gene_name]], log_tranformed=T)
  log_eval_mae = rbind(log_eval_mae,c(lm_mae,rf_mae,svm_mae))
}
# comparing eval_mae and log_mae, svm+log transformed data is better 
# go with svm+log transformed data
estimator = log_svm_estimator 

##################################################
# Network perturbation
# perturb gene expression
get_tree = function(root_name, graph){
  # get the propagation path from roots
  edge_list = get.edgelist(graph)
  edge_list = subset(edge_list,!(edge_list[,2] %in% root_name))
  tree = graph_from_edgelist(edge_list)
  return(tree)
}

dfs_perturb = function(perturbed_gene_name, perturbed_val, tree, estimator,all_node_name,depth=0){
  # heuristic DFS to traverse all leaf nodes
  if(length(perturbed_gene_name)==0|depth==10){
    return(perturbed_val)
  }
  leaf_dist = sapply(perturbed_gene_name,function(x){
    # Depth-first search
    dfs(graph = tree, root =x, dist=T, mode ='out', unreachable=F)$dist
  })
  max_dist = apply(leaf_dist,1,function(x){
    if(all(is.na(x))){
      return(-1)
    }else{
      return(max(x,na.rm=T))
    }
  })
  # only change the 1st level
  next_perturbed_gene = names(which(max_dist==1))
  edgelist = get.edgelist(tree)
  for(n in next_perturbed_gene){
    idx = which(all_node_name==n)
    predictor = estimator[[n]]
    feature_gene = subset(edgelist[,1],edgelist[,2]==n)
    feature_gene_val = perturbed_val[match(feature_gene, all_node_name)]
    feature_gene_val = data.frame(t(feature_gene_val))
    names(feature_gene_val) = feature_gene
    new_n_val = predict(predictor, log(1+feature_gene_val))
    perturbed_val[idx] = exp(new_n_val)-1
  }
  return( dfs_perturb(next_perturbed_gene, perturbed_val, tree, estimator,all_node_name,depth+1) )
}

perturb_grn = function(grn, initial_gene_val, estimator, order=1, target_gene_list = c(), increase_val=0.5){
  candidate_gene = unique(get.edgelist(grn)[,1])
  all_gene = V(grn)$name
  if(length(target_gene_list)==0){
    # Generate All Combinations of n Elements, Taken m at a Time
    all_comb = combn(candidate_gene,order)
  }else{
    all_comb = combn(target_gene_list,order)
  }
  n_comb = dim(all_comb)[2]
  res = list()
  for(i in 1:n_comb){
    print(i)
    curr_node = all_comb[,i]
    curr_idx = match(curr_node, all_gene)
    tree = get_tree(curr_node, grn)
    if(any(!(curr_node %in% V(tree)$name))){
      next
    }
    # case 1: add by 5 (default)
    curr_gene_val = initial_gene_val
    curr_gene_val[curr_idx] = curr_gene_val[curr_idx]+increase_val
    tmp_res = dfs_perturb(curr_node, curr_gene_val, tree, estimator,all_gene)
    curr_res_name = paste0(paste0(curr_node,collapse=','),";add expression")
    names(tmp_res) = all_gene
    res[[curr_res_name]] = tmp_res
    # case 2: mute
    curr_gene_val = initial_gene_val
    curr_gene_val[curr_idx] = 0
    tmp_res = dfs_perturb(curr_node, curr_gene_val, tree, estimator,all_gene)
    curr_res_name = paste0(paste0(curr_node,collapse=','),";mute")
    names(tmp_res) = all_gene
    res[[curr_res_name]] = tmp_res
  }
  return(res)
}

gene_node_name = V(example_graph)$name
current_gene_exp_sub = current_gene_exp[match(gene_node_name, gene_name[,1])]
grn = example_graph
res1 = perturb_grn(grn, current_gene_exp_sub, estimator, order=1)
# long run time warning!!!!
# res2 = perturb_grn(grn, current_gene_exp_sub, estimator, order=2)
# res3 = perturb_grn(grn, current_gene_exp_sub, estimator, order=3)

# example purturbation
get_summary_table <- function(over_exp,mute_exp,purturbed_gene){
  summary_table = as.data.frame(cbind(over_exp,mute_exp))
  summary_table = subset(summary_table,!(rownames(summary_table) %in% purturbed_gene))
  summary_table$logfc = log(summary_table[,1]/summary_table[,2])
  summary_table = summary_table[order(summary_table$logfc,decreasing=T),]
  summary_table = subset(summary_table,summary_table$logfc!=0)
  return(summary_table)
}

eg_res = perturb_grn(grn, current_gene_exp_sub, estimator, order=1
                     , target_gene_list = c("BACH2","RUNX2"))
# ("BACH2","KLF1","RUNX2")

# candidate_gene = unique(get.edgelist(grn)[,1])
# over_exp = eg_res[[1]]
# mute_exp = eg_res[[2]]
over_exp = eg_res[[1]]
mute_exp = eg_res[[2]]
res = get_summary_table(over_exp,mute_exp,c("BACH2","RUNX2"))# ("BACH2","KLF1","RUNX2")







