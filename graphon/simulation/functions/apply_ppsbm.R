

apply_ppsbm = function(case, SEED)
{
  total_time = 50
  
  set.seed(SEED)
  
  if (case==1) network = generate_network1(SEED, total_time)
  if (case==2) network = generate_network2(SEED, total_time)
  if (case==3) network = generate_network3(SEED, total_time)
  if (case==4) network = generate_network4(SEED, total_time)
  
  
  network = del_iso_nodes(network)
  edge_time_mat = network$edge_time_mat
  
  library(ppsbm)
  time.seq = numeric(sum(edge_time_mat<Inf))
  type.seq = numeric(sum(edge_time_mat<Inf))
  current_ind = 1
  for (i in 1:nrow(edge_time_mat)) {
    for (j in 1:ncol(edge_time_mat)) {
      if (edge_time_mat[i,j]<Inf){
        time.seq[current_ind] = edge_time_mat[i,j]
        type.seq[current_ind] = convertNodePair(i, j, n = nrow(edge_time_mat), directed = FALSE)
        current_ind = current_ind+1
      }
    }
  }
  data = list(time.seq=time.seq, type.seq=type.seq, Time=total_time)
  Nijk = statistics(data, nrow(edge_time_mat), K=2^5, directed = FALSE)
  
  # res = mainVEM(data=data, n=nrow(edge_time_mat), Qmin=3, directed=FALSE, method="kernel",
  #               d_part=0, n_perturb=0)
  res = mainVEM(data=list(Nijk=Nijk, Time=total_time), n=nrow(edge_time_mat), 
                Qmin=3, directed=FALSE, method="hist")[[1]]
  res$clusters = mem2clus(apply(res$tau, 2, which.max)) 
  
  return(list(network=network, clus_result=res))
}

