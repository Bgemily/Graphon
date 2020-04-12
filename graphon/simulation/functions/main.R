

main = function(case, SEED, N_clus, N_overclus=N_clus, MaxIter=3, N_subj=1, bw=1)
{
  total_time = 50
  t_vec = seq(0, total_time, 0.05)
  
  if (case==1) network = generate_network1(SEED, total_time)
  if (case==2) network = generate_network2(SEED, total_time)
  if (case==3) network = generate_network3(SEED, total_time)
  if (case==4) network = generate_network4(SEED, total_time)
  
  
  network = del_iso_nodes(network)
  edge_time_mat = network$edge_time_mat
  
  res = do_cluster(edge_time_mat, N_clus, N_overclus, MaxIter, t_vec, bw)
  
  return(list(network=network, clus_result=res))
}



# test
# main(case=2, SEED=67, N_clus=3, N_overclus=4, MaxIter = 1) -> tmp

