get_ARI = function(membership_true, results, truncate=length(membership_true), clusters_old=FALSE)
{
  ARI = vector('numeric', length(results))
  for (i in 1:length(results)) {
    clusters_i = results[[i]]$clusters
    if (clusters_old) clusters_i = results[[i]]$clusters_old

    membership_i = unlist(clusters_i)
    for (j in 1:length(clusters_i)) {
      membership_i[clusters_i[[j]]] = j
    }
    ARI[i] = mclust::adjustedRandIndex(membership_i[1:truncate], membership_true[1:truncate])
  }
  return(ARI)
}